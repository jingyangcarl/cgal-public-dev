// Copyright (c) 2019 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename Conditions>
  class Grouping_segments_2 {
  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;

    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Targets_map = 
      std::map<Size_pair, std::pair<FT, std::size_t> >;
    using Relations_map = 
      std::map<Size_pair, std::pair<int, std::size_t> >;

    Grouping_segments_2() :
    m_tolerance(FT(1) / FT(1000000)),
    m_margin_of_error(FT(0)) 
    { }

    void make_groups(
      const FT max_bound, 
      const std::size_t n, 
      const std::map<std::size_t, Segment_data>& segments,
      const std::vector<FT>& qp_result,
      std::map<FT, Indices>& groups_by_value,
      const Targets_map& targets, 
      const Relations_map& relations = Relations_map()) { 
      
      CGAL_precondition(n > 0);
      CGAL_precondition(max_bound > FT(0));
      CGAL_precondition(qp_result.size() > 0);

      m_conditions.set_margin_of_error(max_bound);
      m_margin_of_error = m_conditions.get_margin_of_error();
      CGAL_assertion(m_margin_of_error > FT(0));
      
      m_groups.clear();
      groups_by_value.clear();
      m_segments_to_groups.clear();
      m_values.clear();

      for (const auto& pair : segments) {
        const auto& seg_data = pair.second;
        const std::size_t seg_index = seg_data.index;
        m_segments_to_groups[seg_index] = -1;
      }

      build_initial_groups(
        n, targets, relations, qp_result);
      build_map_of_values(
        qp_result, segments);

      // Try to assign segments whose orientation has not been optimized 
      // thanks to the regularization process, to an existing group.
      assign_segments_to_groups(segments);
      build_groups_by_value(groups_by_value);
    } 

  private:
    const FT m_tolerance;
    FT m_margin_of_error;
    Conditions m_conditions;
    std::map<std::size_t, int> m_segments_to_groups;
    std::map<std::size_t, Indices> m_groups;
    std::map<int, FT> m_values;
    
    void build_initial_groups(
      const std::size_t n,
      const Targets_map& targets, 
      const Relations_map& relations,
      const std::vector<FT>& qp_result) {
      
      std::size_t g = 0;
      auto rel_it = relations.begin();

      for (const auto& target : targets) {
        const std::size_t i = target.first.first;
        const std::size_t j = target.first.second;
        const std::size_t p = target.second.second;

        int r = 0;
        if (rel_it != relations.end()) {
          CGAL_assertion(rel_it->second.second == p);
          r = rel_it->second.first;
        }
        CGAL_assertion(r == 0 || r == 1);

        if (CGAL::abs(qp_result[n + p]) >= m_tolerance) { 
          if (rel_it != relations.end()) ++rel_it;
          continue;
        }

        const int g_i = m_segments_to_groups[i];
        const int g_j = m_segments_to_groups[j];
        const int status = check_group_status(g_i, g_j);

        switch (status) {
          case -1: 
            break;
          case 1: {
            if (r == 0) create_single_group(i, j, g);
            else create_separate_groups(i, j, g); 
            break;
          }
          case 2: {
            if (r == 0) assign_segment_to_group(i, j);
            else create_new_group(i, g);
            break; 
          }
          case 3: {
            if (r == 0) assign_segment_to_group(j, i);
            else create_new_group(j, g); 
            break; 
          }
          case 4: {
            if (r == 0) merge_two_groups(g_i, g_j); 
            break; 
          }
        }
        if (rel_it != relations.end()) 
          ++rel_it;
      }
    }

    void build_map_of_values(
      const std::vector<FT>& qp_result,
      const std::map<std::size_t, Segment_data>& segments) {
      
      for (const auto& segment_to_group : m_segments_to_groups) {
        int g_i = segment_to_group.second;
        if (g_i != -1 && (m_values.find(g_i) == m_values.end())) {
          const std::size_t seg_index = segment_to_group.first;

          const auto& seg_data = segments.at(seg_index);
          const FT value = m_conditions.reference(
            seg_data, qp_result[seg_index]);
          
          // Check if the angle that seems to be associated to this group 
          // of segments is not too close to another value.
          int g_j = -1;
          for (const auto& pair : m_values) {
            if (CGAL::abs(pair.second - value) < m_margin_of_error) 
              g_j = pair.first;
          }

          if (g_j == -1) m_values[g_i] = value;
          else merge_two_groups(g_j, g_i);
        }
      }
    }

    void assign_segments_to_groups(
      const std::map<std::size_t, Segment_data>& segments) {
      
      for (const auto& segment_to_group : m_segments_to_groups) {
        int g_i = segment_to_group.second;
        if (g_i == -1) {
          const std::size_t seg_index = segment_to_group.first;

          const auto& seg_data = segments.at(seg_index);
          const FT value = m_conditions.reference(seg_data, 0);
          
          int g_j = -1;
          for (const auto& pair : m_values) {
            const FT value_j = pair.second;
            const int g_index = pair.first;

            g_j = m_conditions.group_index(
              value, value_j, g_index);
            if (g_j != -1) break;
          }

          if (g_j == -1) { 
            m_values.size() > 0 ? g_i = m_values.rbegin()->first + 1 : g_i = 0;
            m_values[g_i] = value;
          } else g_i = g_j;

          m_segments_to_groups[seg_index] = g_i;
          m_groups[g_i].push_back(seg_index);
        }
      }
    }

    void build_groups_by_value(
      std::map<FT, Indices>& groups_by_value) {
      for (const auto& pair : m_values) {
        const FT value = pair.second;
        if (groups_by_value.find(value) == groups_by_value.end()) 
          groups_by_value[value] = Indices();
      }

      for (const auto& segment_to_group : m_segments_to_groups) {
        const FT value = m_values.at(segment_to_group.second);     
        if (groups_by_value.find(value) != groups_by_value.end()) 
          groups_by_value[value].push_back(segment_to_group.first);
      }
    } 

    int check_group_status(const int g_i, const int g_j) const {
      if (g_i == -1 && g_j == -1) return 1;
      if (g_i == -1 && g_j != -1) return 2;
      if (g_i != -1 && g_j == -1) return 3;
      if (g_i != -1 && g_j != -1 && g_i != g_j) return 4;
      return -1;
    }

    void create_single_group(
      const std::size_t i, const std::size_t j, 
      std::size_t& g) {

      m_segments_to_groups[i] = g;
      m_segments_to_groups[j] = g;
      m_groups[g].push_back(i);
      m_groups[g].push_back(j); 
      ++g;
    }

    void create_separate_groups(
      const std::size_t i, const std::size_t j, 
      std::size_t& g) {
      
      create_new_group(i, g);
      create_new_group(j, g);
    }

    void assign_segment_to_group(
      const std::size_t i, const std::size_t j) {
      
      const int g_j = m_segments_to_groups[j];
      m_segments_to_groups[i] = g_j;
      m_groups[g_j].push_back(i);
    }

    void create_new_group(
      const std::size_t i, 
      std::size_t& g) {
      
      m_segments_to_groups[i] = g;
      m_groups[g].push_back(i);
      ++g;
    }

    void merge_two_groups(
      const int g_i, const int g_j) {
      
      for (const std::size_t index : m_groups[g_j]) {
        m_segments_to_groups[index] = g_i;
        m_groups[g_i].push_back(index);
      }
      m_groups[g_j].clear();
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_GROUPING_SEGMENTS_2_H
