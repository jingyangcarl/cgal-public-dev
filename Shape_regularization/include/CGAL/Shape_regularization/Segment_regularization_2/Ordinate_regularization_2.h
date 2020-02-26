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

#ifndef CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>
#include <CGAL/Shape_regularization/internal/Ordinate_conditions_2.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularization2DReg

    \brief An ordinate-based regularization type on a set of 2D segments that preserves 
    collinearity relationship.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`.

    \cgalModels `RegularizationType`
  */
  template<
  typename GeomTraits, 
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Ordinate_regularization_2 {
  public:

    /// \name Types
    /// @{
    
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// \cond SKIP_IN_MANUAL
    using Point_2 = typename Traits::Point_2;
    using Vector_2  = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;

    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Conditions = typename internal::Ordinate_conditions_2<Traits>;
    using Grouping = internal::Grouping_segments_2<Traits, Conditions>;

    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;

    using Targets_map = 
      std::map<Size_pair, std::pair<FT, std::size_t> >;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param d_max
      max distance value in meters

      \param segment_map
      an instance of `SegmentMap` that maps an item 
      from `input_range` to `GeomTraits::Segment_2`

      \pre `input_range.size() > 1`
      \pre `d_max >= 0`
    */
    Ordinate_regularization_2 (
      InputRange& input_range,
      const FT d_max = FT(1) / FT(10),
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_d_max(CGAL::abs(d_max)),
    m_segment_map(segment_map),
    m_num_modified_segments(0) { 

      CGAL_precondition(input_range.size() > 1);
      CGAL_precondition(d_max >= FT(0));

      if (d_max < FT(0)) {
        std::cout << 
          "WARNING: The max ordinate bound has to be within [0, +inf)! Setting to 0." 
        << std::endl;
        m_d_max = FT(0);
      }
    }
    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief implements `RegularizationType::target_value()`.

      This function calculates the target value between 2 segments, which are
      direct neighbors to each other.

      \param query_index_i
      index of the first segment

      \param query_index_j
      index of the second segment

      \pre `query_index_i >= 0 && query_index_i < input_range.size()`
      \pre `query_index_j >= 0 && query_index_j < input_range.size()`
    */

    FT target_value(
      const std::size_t query_index_i, 
      const std::size_t query_index_j) {

      if( m_segments.size() == 0) 
        return FT(0);

      CGAL_precondition(m_segments.size() > 1);
      CGAL_precondition(m_segments.find(query_index_i) != m_segments.end());
      CGAL_precondition(m_segments.find(query_index_j) != m_segments.end());

      const std::size_t i = query_index_i;
      const std::size_t j = query_index_j;

      const auto& s_i = m_segments.at(i);
      const auto& s_j = m_segments.at(j);

      const FT tar_val = 
        s_i.ref_coords.y() - s_j.ref_coords.y();
      if (CGAL::abs(tar_val) < bound(i) + bound(j))
        m_targets[std::make_pair(i, j)] = tar_val;
      return tar_val;
    }

    /*!
      \brief implements `RegularizationType::bound()`.

      This function returns `theta_max`.
    */
    FT bound(const std::size_t) const {
      return m_d_max;
    }

    /*!
      \brief implements `RegularizationType::update()`.

      This function applies new positions computed by the QP solver 
      to the initial segments.

      \param result
      a vector with updated segment positions.

      \pre `result.size() > 0`
    */
    void update(const std::vector<FT>& result) {
      CGAL_precondition(result.size() > 0);

      Targets_map targets;
      std::map<FT, Indices> collinear_groups;
      std::map<std::size_t, Segment_data> segments;

      CGAL_precondition(m_targets.size() > 0);
      for (const auto& group : m_groups) {
        if (group.size() < 2) continue; 

        targets.clear(); segments.clear();
        build_grouping_data(group, segments, targets);

        collinear_groups.clear();
        if (segments.size() > 0) {
          const std::size_t n = m_input_range.size();

          m_grouping.make_groups(
            m_d_max, n, segments, result, 
            collinear_groups, targets);
          translate_collinear_segments(collinear_groups);
        }
      }
    }
    /// @}

    /// \name Utilities
    /// @{ 

    /*!
      \brief inserts a group of segments from `input_range`.

      \tparam ItemRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

      \tparam IndexMap 
      must be an `LvaluePropertyMap` whose key type is the value type of `ItemRange`
      and value type is `std::size_t`.

      \param item_range
      an instance of ItemRange

      \param index_map
      an instance of IndexMap that returns an index stored in the `item_range` 
      of the segment in the `input_range`

      \pre `item_range.size() > 1`
    */
    template<
    typename ItemRange, 
    typename IndexMap = CGAL::Identity_property_map<std::size_t> >
  	void add_group(
      const ItemRange& item_range, 
      const IndexMap index_map = IndexMap()) { 
      
      CGAL_precondition(item_range.size() > 1);
      if (item_range.size() < 2) return;
      
      Indices group;
      group.reserve(item_range.size());
      for (const auto& item : item_range) {
        const std::size_t seg_index = get(index_map, item);
        group.push_back(seg_index);
      }
      
      m_groups.push_back(group);
      update_segment_data(group);
    }

    /*!
      \brief returns number of modifed segments`.
    */
    std::size_t number_of_modified_segments() const {
      return m_num_modified_segments;
    }

    /// @}

  private:
    Input_range& m_input_range;
    FT m_d_max;
    const Segment_map m_segment_map;
    std::map<std::size_t, Segment_data> m_segments;
    std::map<Size_pair, FT> m_targets;
    Grouping m_grouping;
    std::vector<Indices> m_groups;
    std::size_t m_num_modified_segments;

    void update_segment_data(
      const Indices& paral_gr) {
      if (paral_gr.size() < 2) return;

      Point_2 frame_origin;
      for(std::size_t i = 0; i < paral_gr.size(); ++i) {
        const std::size_t seg_index = paral_gr[i];

        CGAL_precondition(m_segments.find(seg_index) == m_segments.end());
        if(m_segments.find(seg_index) != m_segments.end())
          continue;

        const auto& seg = get(m_segment_map, *(m_input_range.begin() + seg_index));
        Segment_data seg_data(seg, seg_index);

        if (i == 0)
          frame_origin = seg_data.barycenter;

        seg_data.ref_coords = internal::transform_coordinates_2(
                seg_data.barycenter, frame_origin, seg_data.orientation);
        m_segments.emplace(seg_index, seg_data);
      } 
    }

    void build_grouping_data(
      const std::vector <std::size_t> & group,
      std::map <std::size_t, Segment_data> & segments,
      Targets_map & targets) {
      
      for (const std::size_t it : group) {
        const std::size_t seg_index = it;

        CGAL_precondition(m_segments.find(seg_index) != m_segments.end());
        const Segment_data& seg_data = m_segments.at(seg_index);

        segments.emplace(seg_index, seg_data);
        std::size_t tar_index = 0;

        for(const auto & ti : m_targets) {
          const std::size_t seg_index_tar_i = ti.first.first;
          const std::size_t seg_index_tar_j = ti.first.second;
          const FT tar_val = ti.second;

          if (seg_index_tar_i == seg_index) {
            targets[std::make_pair(seg_index_tar_i, seg_index_tar_j)] = std::make_pair(tar_val, tar_index);
          }

          ++tar_index;
        }
      }
    }

    void translate_collinear_segments(
      const std::map <FT, std::vector<std::size_t>> & collinear_groups_by_ordinates) {
      for (const auto & mi : collinear_groups_by_ordinates) {
        const FT dt = mi.first;
        const std::vector<std::size_t> & group = mi.second;
        int l_index = find_longest_segment(group);
        CGAL_postcondition(l_index >= 0);
        if(l_index < 0) {
          std::cerr << "Cannot translate collinear segments! Cannot find the longest segment!" << std::endl;
          return;
        }

        CGAL_precondition(m_segments.find(l_index) != m_segments.end());
        const auto& l_data = m_segments.at(l_index);

        FT new_difference = dt - l_data.ref_coords.y();
        set_difference(l_index, new_difference);

        const FT l_a = l_data.a;
        const FT l_b = l_data.b;
        const FT l_c = l_data.c;
        const auto & l_direction = l_data.direction;

        // Translate the other segments, so that they rest upon the line ax + by + c = 0.
        for (const std::size_t it : group) {
          if (it != static_cast<std::size_t>(l_index)) {
            CGAL_precondition(m_segments.find(it) != m_segments.end());
            const Segment_data & seg_data = m_segments.at(it);

            new_difference = dt - seg_data.ref_coords.y();
            set_difference(it, new_difference, l_a, l_b, l_c, l_direction);
          }
        }
      }
    }

    int find_longest_segment(
      const std::vector<std::size_t> & group) const {
      FT l_max = -FT(1000000000000);
      int l_index = -1;

      for (const std::size_t it : group) {
        const FT seg_length = m_segments.at(it).length;

        if (l_max < seg_length) {
          l_max = seg_length;
          l_index = it;
        }
      }

      return l_index;
    }

    void set_difference(
      const int i, const FT new_difference) {
      const FT difference = new_difference;
      Segment_data & seg_data = m_segments.at(i);

      const auto & direction = seg_data.direction;
      const Vector_2 final_normal = Vector_2(-direction.y(), direction.x());

      const auto &source = seg_data.segment.source();
      const auto &target = seg_data.segment.target();

      Point_2 new_source = Point_2(source.x() + difference * final_normal.x(), source.y() + difference * final_normal.y());
      Point_2 new_target = Point_2(target.x() + difference * final_normal.x(), target.y() + difference * final_normal.y());
      
      const FT bx = (new_source.x() + new_target.x()) / FT(2);
      const FT by = (new_source.y() + new_target.y()) / FT(2);

      m_input_range[i] = Segment_2(new_source, new_target);
      seg_data.c = -seg_data.a * bx - seg_data.b * by;

      ++m_num_modified_segments;
    }

    void set_difference(
      const int i, const FT new_difference, const FT a, const FT b, const FT c, const Vector_2 &direction) {
      FT difference = new_difference;
      auto & seg_data = m_segments.at(i);

      seg_data.direction = direction;
      if (seg_data.direction.y() < FT(0) || (seg_data.direction.y() == FT(0) && seg_data.direction.x() < FT(0))) 
        seg_data.direction = -seg_data.direction;

      Vector_2 final_normal = Vector_2(-seg_data.direction.y(), seg_data.direction.x());
      FT x1, x2, y1, y2;

      const auto &source = seg_data.segment.source();
      const auto &target = seg_data.segment.target();

      if (CGAL::abs(seg_data.direction.x()) > CGAL::abs(seg_data.direction.y())) {
        x1 = source.x() + difference * final_normal.x();
        x2 = target.x() + difference * final_normal.x(); 

        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;
      } 
      else {    
        y1 = source.y() + difference * final_normal.y();
        y2 = target.y() + difference * final_normal.y();

        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }

      const Point_2 new_source = Point_2(x1, y1);
      const Point_2 new_target = Point_2(x2, y2);
      m_input_range[i] = Segment_2(new_source, new_target);

      ++m_num_modified_segments;
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ORDINATE_REGULARIZATION_2_H
