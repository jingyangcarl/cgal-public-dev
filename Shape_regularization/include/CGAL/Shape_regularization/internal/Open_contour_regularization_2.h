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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Andreas Fabri
//

#ifndef CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Shape_regularization/enum.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Open_contour_regularization_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    struct Segment_wrapper_2 {

      Segment_2 segment;
      std::size_t index = std::size_t(-1);
      std::size_t object_index = std::size_t(-1);
      bool is_valid_direction = false;
      bool is_used = false;
    };

    using FT_pair = std::pair<FT, FT>;
    using Segments_2 = std::vector<Segment_2>;
    using Segment_wrappers_2 = std::vector<Segment_wrapper_2>;

    Open_contour_regularization_2(
      Input_range& input_range,
      Point_map point_map,
      const FT angle_threshold_2 = FT(5),
      const bool verbose = true) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_angle_threshold_2(angle_threshold_2),
    m_verbose(verbose) { 
      
      CGAL_precondition(input_range.size() >= 2);
    }

    void set_principal_directions(
      const std::vector<std::size_t>& directions) {

      CGAL_precondition(
        directions.size() == m_input_range.size() - 1);
      m_group = directions;
    }

    void estimate_principal_directions(
      const Direction_type direction_type,
      const FT min_length_2,
      const FT max_angle_2) {

      switch (direction_type) {
        case Direction_type::LONGEST:
          set_longest_direction();
          break;
        case Direction_type::MULTIPLE:
          set_multiple_directions(
            min_length_2, max_angle_2);
          break;
        default:
          set_longest_direction();
          break;
      };
    }

    template<typename OutputIterator>
    void regularize(
      OutputIterator contour,
      const FT max_ordinate_2) {
      
      Examples::Saver<Traits> saver;

      bool success = false;
      rotate_contour();
      if (m_verbose)
        saver.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/rotated");

      exit(EXIT_SUCCESS);

      /*
      success = optimize_contour(max_ordinate_2);
      if (!success) return;
      if (m_verbose)
        saver.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/optimized");

      success = connect_contour(max_ordinate_2);
      if (!success) return;
      if (m_verbose)
        saver.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/connected");

      update_input(contour); */
    }

    std::size_t number_of_principal_directions() const {

      std::set<std::size_t> unique;
      for (std::size_t index : m_group)
        if (index != std::size_t(-1))
          unique.insert(index);
      return unique.size();
    }

  private:
    Input_range& m_input_range;
    Point_map m_point_map;

    const FT m_angle_threshold_2;
    const bool m_verbose;

    std::vector<FT_pair> m_bounds;
    std::vector<Segment_2> m_longest;
    std::vector<std::size_t> m_group;

    std::vector<Segment_wrapper_2> m_wraps;

    void set_longest_direction() {

      m_bounds.clear(); m_bounds.resize(1);
      m_bounds[0] = std::make_pair(FT(45), FT(45));

      m_longest.clear(); m_longest.resize(1);
      m_longest[0] = compute_longest_segment();

      m_group.clear();
      m_group.resize(m_input_range.size(), 0); // 0 is the index of the direction in the m_longest
    }

    Segment_2 compute_longest_segment() const {

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      const std::size_t n = m_input_range.size();
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        const FT length = CGAL::squared_distance(source, target);
        if (length > max_length) {

          longest = i;
          max_length = length;
        }
      }

      const std::size_t i = longest;
      const std::size_t ip = i + 1;
      const auto& source = get(m_point_map, *(m_input_range.begin() + i));
      const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

      return Segment_2(source, target);
    }

    void set_multiple_directions(
      const FT min_length_2,
      const FT max_angle_2) {
      
      create_segments_from_input(min_length_2);
      estimate_initial_directions(max_angle_2);

      if (m_longest.size() == 0) {
        set_longest_direction();
      } else {
        unify_along_contours();
        correct_directions();
        readjust_directions();
      }

      if (m_verbose) {
        std::cout << "* groups: ";
        for (std::size_t group_index : m_group)
          std::cout << group_index << " ";
        std::cout << std::endl;
      }
    }

    void create_segments_from_input(
      const FT min_length_2) {

      const std::size_t n = m_input_range.size();
      m_wraps.clear();
      m_wraps.reserve(n);
      
      Segment_wrapper_2 wrap;
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;
        
        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        wrap.index = i;
        wrap.segment = Segment_2(source, target);
        wrap.is_valid_direction = 
          is_valid_principal_direction(wrap.segment, min_length_2);
        m_wraps.push_back(wrap);
      }
      CGAL_assertion(m_wraps.size() == m_input_range.size() - 1);
    }

    bool is_valid_principal_direction(
      const Segment_2& segment,
      const FT min_length_2) const {
      
      CGAL_precondition(min_length_2 >= FT(0));
      return internal::length_2(segment) >= min_length_2 * FT(2);
    }

    void estimate_initial_directions(
      const FT max_angle_2) {

      std::vector<std::size_t> longest_to_short;
      sort_segments_by_length(m_wraps, longest_to_short);

      m_bounds.clear(); m_longest.clear(); m_group.clear();
      m_group.resize(longest_to_short.size(), std::size_t(-1));

      std::size_t group_index = 0;
      std::size_t query_index = std::size_t(-1);
      do {
        query_index = find_next_longest_segment(longest_to_short);
        if (query_index != std::size_t(-1))
          set_next_longest_direction(
            max_angle_2, query_index, group_index);
        ++group_index;
      } while (query_index != std::size_t(-1));
    }

    void sort_segments_by_length(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& sorted) const {

      sorted.clear();
      sorted.reserve(wraps.size());
      for (std::size_t i = 0; i < wraps.size(); ++i)
        sorted.push_back(i);

      std::sort(sorted.begin(), sorted.end(), 
      [&wraps](const std::size_t i, const std::size_t j) -> bool { 
        
        const FT length_1 = wraps[i].segment.squared_length();
        const FT length_2 = wraps[j].segment.squared_length();
        return length_1 > length_2;
      });
    }

    std::size_t find_next_longest_segment(
      const std::vector<std::size_t>& longest_to_short) const {

      std::size_t longest = std::size_t(-1);
      for (std::size_t i = 0; i < longest_to_short.size(); ++i) {
        const std::size_t wrap_index = longest_to_short[i];
        const auto& wrap = m_wraps[wrap_index];

        if (is_valid_wrap(wrap)) {
          longest = wrap_index; break;
        }
      }
      return longest;
    }

    bool is_valid_wrap(
      const Segment_wrapper_2& wrap) const {  
      return !wrap.is_used && wrap.is_valid_direction;
    }

    void set_next_longest_direction(
      const FT max_angle_2,
      const std::size_t query_index,
      const std::size_t group_index) {

      CGAL_assertion(query_index != std::size_t(-1));
      CGAL_assertion(group_index != std::size_t(-1));
      
      // Set current longest direction.
      auto& longest = m_wraps[query_index];
      m_group[query_index] = group_index;
      longest.is_used = true;

      for (auto& wrap : m_wraps) {
        if (wrap.index == query_index) // skip longest
          continue;

        // Check if another wrap satisifes the conditions.
        if (is_valid_wrap(wrap)) { 
          if (does_satisify_angle_conditions(
            max_angle_2, longest.segment, wrap.segment)) {
            
            m_group[wrap.index] = group_index;
            wrap.is_used = true;
          }
        }
      }

      // Set internals.
      m_longest.push_back(longest.segment);
      m_bounds.push_back(std::make_pair(FT(45), FT(45)));
    }

    bool does_satisify_angle_conditions(
      const FT max_angle_2,
      const Segment_2& longest,
      const Segment_2& segment) const {

      CGAL_precondition(
        max_angle_2 >= FT(0) && max_angle_2 <= FT(90));
      const FT bound_min = max_angle_2;
      const FT bound_max = FT(90) - bound_min;

      const FT angle_2 = CGAL::abs(
        internal::angle_2_degrees(longest, segment));
      return (angle_2 <= bound_min) || (angle_2 >= bound_max);
    }

    void unify_along_contours() {

      const std::size_t n = m_wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrap = m_wraps[i];
        if (wrap.is_used) continue;
        
        std::size_t im = std::size_t(-1);
        if (i > 0) im = i - 1;
        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;

        bool stop = false;
        std::size_t max_count = 0;
        do {

          if (im != std::size_t(-1) && m_wraps[im].is_used) {
            m_group[i] = m_group[im];
            m_wraps[i].is_used = true;
            break;
          }

          if (ip != std::size_t(-1) && m_wraps[ip].is_used) {
            m_group[i] = m_group[ip]; 
            m_wraps[i].is_used = true;
            break;
          }

          if (im != std::size_t(-1) && im > 0) 
            im = im - 1;
          if (ip != std::size_t(-1) && ip < n - 1) 
            ip = ip + 1;

          if (im == 0 || ip == n - 1) // fix this, I skip 0 and last
            stop = true;
          ++max_count;

        } while (!stop && max_count < n);
        if (stop || max_count >= n) {
          std::cerr << 
            "Warning: revert back to the first direction!" << std::endl;
          m_group[i] = 0;
        }
      }
    }

    void correct_directions() {

      const std::size_t n = m_wraps.size();
      std::vector<std::size_t> clean;
      clean.reserve(n);
      
      for (std::size_t i = 0; i < n; ++i) {
        
        std::size_t im = std::size_t(-1); // fix this, I skip 0 and last
        if (i > 0) im = i - 1;
        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;

        std::size_t gm = std::size_t(-1);
        if (im != std::size_t(-1)) gm = m_group[im];
        std::size_t gi = std::size_t(-1);
        if (i != std::size_t(-1)) gi = m_group[i];
        std::size_t gp = std::size_t(-1);
        if (ip != std::size_t(-1)) gp = m_group[ip];

        if (gm != std::size_t(-1) && gm == gp && gi != gm)
          clean.push_back(gm);
        else
          clean.push_back(gi);
      }
      m_group = clean;
    }

    void readjust_directions() {

      std::vector<FT> angles, counts;
      create_average_angles(angles, counts);

      CGAL_assertion(angles.size() == counts.size());
      CGAL_assertion(angles.size() == m_longest.size());

      for (std::size_t k = 0; k < angles.size(); ++k) {
        CGAL_assertion(counts[k] != FT(0));
        angles[k] /= counts[k];

        const FT angle_deg = angles[k];
        rotate(angle_deg, FT(0), m_longest[k]);
      }
    }

    void create_average_angles(
      std::vector<FT>& angles,
      std::vector<FT>& counts) const {

      CGAL_assertion(m_longest.size() > 0);

      angles.clear();
      angles.resize(m_longest.size(), FT(0));

      counts.clear();
      counts.resize(m_longest.size(), FT(0));

      for (std::size_t i = 0; i < m_wraps.size(); ++i) {
        const auto& wrap = m_wraps[i];
        if (!wrap.is_valid_direction) continue;

        const std::size_t group_index = m_group[i];
        CGAL_assertion(group_index != std::size_t(-1));

        const auto& si = m_longest[group_index];
        const auto& sj = wrap.segment;

        const auto di = internal::direction_2(si);
        const auto dj = internal::direction_2(sj);

        const FT oi = internal::orientation_2(di);
        const FT oj = internal::orientation_2(dj);

        const FT mes_ij = oi - oj;
        const double mes90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

        const FT to_lower = FT(90) *  static_cast<FT>(mes90)          - mes_ij;
        const FT to_upper = FT(90) * (static_cast<FT>(mes90) + FT(1)) - mes_ij;

        const FT angle = 
          CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

        angles[group_index] += angle;
        counts[group_index] += FT(1);
      }
    }

    void rotate(
      const FT angle_2, // in degrees
      const FT ref_angle_2, // in degrees
      Segment_2& segment) const {

      FT angle = angle_2;
      if (angle < FT(0)) angle = angle + ref_angle_2;
      else if (angle > FT(0)) angle = angle - ref_angle_2;

      Point_2 source = segment.source();
      Point_2 target = segment.target();

      const Point_2 b = internal::middle_point_2(source, target);
      const FT angle_rad = angle * static_cast<FT>(CGAL_PI) / FT(180);

      internal::rotate_point_2(angle_rad, b, source);
      internal::rotate_point_2(angle_rad, b, target);

      segment = Segment_2(source, target);
    }

    void rotate_contour() {

      CGAL_assertion(m_group.size() == m_wraps.size());
      for (std::size_t i = 0; i < m_wraps.size(); ++i) {
        const std::size_t group_index = m_group[i];
        if (group_index == std::size_t(-1))
          continue;

        auto& wrap = m_wraps[i];
        const auto& longest = m_longest[group_index];
        const auto& bounds = m_bounds[group_index];
        rotate_segment(
          bounds, longest, wrap.segment);
      }
    }

    void rotate_segment(
      const FT_pair& bounds,
      const Segment_2& longest, 
      Segment_2& segment) {

      const FT angle_deg = internal::compute_angle_2(longest, segment);
      const FT converted = CGAL::abs(convert_angle_2(angle_deg));

      if (converted <= bounds.first)
        rotate(angle_deg, FT(180), segment); // parallel case
      if (converted >= bounds.second)
        rotate(angle_deg, FT(90), segment); // orthogonal case
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H
