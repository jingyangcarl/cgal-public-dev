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
      
      if (m_wraps.size() < 2) return;
      Examples::Saver<Traits> saver;

      bool success = false;
      rotate_contour();
      if (m_verbose)
        saver.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/rotated");

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

      update_input(contour);
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

    bool optimize_contour(
      const FT max_ordinate_2) {

      // Clean.
      remove_zero_length_segments(m_wraps);
      if (m_wraps.size() == 0) return false;

      // Create groups of collinear segments.
      std::vector<Segment_wrappers_2> groups;
      create_consecutive_groups(m_wraps, groups);
      if (m_verbose)
        std::cout << 
          "* number of consecutive groups = " << groups.size() << std::endl;

      // Optimize all groups with at least two collinear segments.
      std::size_t count = 0;
      for (auto& group : groups) {
        if (group.size() > 1) {
          optimize_group(max_ordinate_2, group);
          ++count;
        }
      }
      if (m_verbose)
        std::cout << 
          "* number of optimized groups = " << count << std::endl;

      // Update segments.
      create_segments_from_groups(groups, m_wraps);
      if (m_wraps.size() == 0) return false;
      return true;
    }

    void remove_zero_length_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<Segment_wrapper_2> clean;
      for (const auto& wrap : wraps)
        if (wrap.segment.squared_length() > internal::tolerance<FT>())
          clean.push_back(wrap);
      wraps = clean;
    }

    void create_consecutive_groups(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<Segment_wrappers_2>& groups) const {

      groups.clear();
      for (auto& wrap : wraps)
        wrap.is_used = false;
      std::vector<Segment_wrapper_2> group;

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n - 1; ++i) {
        auto& wrapi = wraps[i];
        if (wrapi.is_used) continue;
        
        group.clear(); 
        wrapi.is_used = true;
        group.push_back(wrapi);
        
        const std::size_t ip = i + 1;
        for (std::size_t j = ip; j < n; ++j) {
          auto& wrapj = wraps[j];
          const FT angle_2 = CGAL::abs(
            internal::angle_2_degrees(wrapi.segment, wrapj.segment));

          if (angle_2 <= m_angle_threshold_2) {
            wrapj.is_used = true;
            group.push_back(wrapj); 
          } else break;
        }
        groups.push_back(group);
      }

      // Handle the last segment in case it is not included in any group.
      auto& last = wraps[n - 1];
      if (!last.is_used) {
        last.is_used = true;

        group.clear();
        group.push_back(last);
        groups.push_back(group);
      }
    }

    void optimize_group(
      const FT max_ordinate_2,
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<std::size_t> longest_to_short;
      sort_segments_by_length(wraps, longest_to_short);

      std::vector<Segment_wrappers_2> groups;
      create_collinear_groups(
        max_ordinate_2, wraps, longest_to_short, groups);

      /*
      if (m_verbose)
        std::cout << 
          "* number of collinear groups = " << groups.size() << std::endl; */

      const std::size_t before = wraps.size();
      std::vector<Line_2> lines;
      create_lines(groups, lines);
      move_segments_toward_lines(lines, wraps);
      update_group(wraps);
      const std::size_t after = wraps.size();
      CGAL_assertion(after >= before);

      /*
      if (m_verbose)
        std::cout << 
          "* segments before/after = " << before << "/" << after << std::endl; */
    }

    void create_collinear_groups(
      const FT max_ordinate_2,
      std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& seeds,
      std::vector<Segment_wrappers_2>& groups) const {

      groups.clear();
      for (auto& wrap : wraps)
        wrap.is_used = false;
      std::vector<Segment_wrapper_2> group;

      std::size_t group_index = 0;
      const int n = static_cast<int>(wraps.size());
      for (const std::size_t seed : seeds) {
        const int i = static_cast<int>(seed);

        auto& wrapi = wraps[i];
        if (wrapi.is_used) continue;
        
        group.clear(); 
        wrapi.is_used = true;
        wrapi.object_index = group_index;
        group.push_back(wrapi);

        const auto source = internal::middle_point_2(
          wrapi.segment.source(), wrapi.segment.target());

        // Traverse forward.
        const int ip = i + 1;
        if (i < n - 1 && !wraps[ip].is_used) {
          int j = ip;
          while (j < n) {
            auto& wrapj = wraps[j];
            if (wrapj.is_used) break;

            if (does_satisfy_ordinate_conditions(
              max_ordinate_2, source, wrapj.segment)) {

              wrapj.is_used = true;
              wrapj.object_index = group_index;
              group.push_back(wrapj);
            } else break;
            ++j;
          }
        }

        // Traverse backward.
        const int im = i - 1;
        if (i > 0 && !wraps[im].is_used) {
          int j = im;
          while (j >= 0) {
            auto& wrapj = wraps[j];
            if (wrapj.is_used) break;
        
            if (does_satisfy_ordinate_conditions(
              max_ordinate_2, source, wrapj.segment)) {

              wrapj.is_used = true;
              wrapj.object_index = group_index;
              group.push_back(wrapj);
            } else break;
            --j;
          }
        }

        groups.push_back(group);
        ++group_index;
      }
    }

    bool does_satisfy_ordinate_conditions(
      const FT max_ordinate_2,
      const Point_2& source,
      const Segment_2& segment) const {

      CGAL_assertion(max_ordinate_2 >= FT(0));
      const Line_2 line = Line_2(
        segment.source(), segment.target());
      const auto target = line.projection(source);
      const Segment_2 proj = Segment_2(source, target);
      const FT distance = internal::length_2(proj);
      return distance <= max_ordinate_2;
    }

    void create_lines(
      const std::vector<Segment_wrappers_2>& groups,
      std::vector<Line_2>& lines) const {

      CGAL_assertion(groups.size() > 0);
      lines.clear();
      lines.reserve(groups.size());

      for (const auto& group : groups) {
        const Segment_2 segment = find_weighted_segment(group);
        const Line_2 line = Line_2(segment.source(), segment.target());
        lines.push_back(line);
      }
    }

    Segment_2 find_weighted_segment(
      const std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<FT> weights;
        compute_distance_weights(wraps, weights);
      const Segment_2 ref_segment = 
        find_central_segment(wraps);
      const Segment_2 weighted = 
        compute_weighted_segment(wraps, weights, ref_segment);
      if (weighted.source() == weighted.target())
        return ref_segment;
      return weighted;
    }

    void compute_distance_weights(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT>& weights) const {

      CGAL_assertion(wraps.size() > 0);
      weights.clear();
      weights.reserve(wraps.size());

      FT sum_distance = FT(0);
      for (const auto& wrap : wraps) {
        const FT distance = internal::length_2(wrap.segment);
        sum_distance += distance;
        weights.push_back(distance);
      }

      CGAL_assertion(sum_distance > FT(0));
      for (auto& weight : weights)
        weight /= sum_distance;
      CGAL_assertion(weights.size() == wraps.size());
    }

    Segment_2 find_central_segment(
      const std::vector<Segment_wrapper_2>& wraps) const {

      Point_2 source, target;
      FT x1 = FT(0), y1 = FT(0);
      FT x2 = FT(0), y2 = FT(0);
      for (const auto& wrap : wraps) {
        x1 += wrap.segment.source().x();
        x2 += wrap.segment.target().x();

        y1 += wrap.segment.source().y();
        y2 += wrap.segment.target().y();
      }

      CGAL_assertion(wraps.size() > 0);
      const FT size = static_cast<FT>(wraps.size());
      x1 /= size; y1 /= size;
      x2 /= size; y2 /= size;

      source = Point_2(x1, y1);
      target = Point_2(x2, y2);

      if (source == target)
        return find_longest_segment(wraps);
      return Segment_2(source, target);
    }

    Segment_2 find_longest_segment(
      const std::vector<Segment_wrapper_2>& wraps) const {

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        const FT length = wrap.segment.squared_length();
        if (length > max_length) {
          longest = i;
          max_length = length;
        }
      }
      return wraps[longest].segment;
    }

    Segment_2 compute_weighted_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<FT>& weights,
      const Segment_2& ref_segment) const {

      const auto& sref = ref_segment.source();
      const auto& tref = ref_segment.target();

      const auto center = 
        internal::middle_point_2(sref, tref);

      CGAL_assertion(weights.size() == wraps.size());
      Vector_2 dir = Vector_2(FT(0), FT(0));
      for (std::size_t i = 0; i < weights.size(); ++i) {  
        const FT weight = weights[i];

        const auto& wrap = wraps[i];
        const Line_2 line = Line_2(
          wrap.segment.source(), wrap.segment.target());
        const Point_2 proj = line.projection(center);

        const Vector_2 v = Vector_2(center, proj);
        dir += v * weight;
      }

      const Point_2 source = sref + dir;
      const Point_2 target = tref + dir;
      return Segment_2(source, target);
    }

    void move_segments_toward_lines(
      const std::vector<Line_2>& lines,
      std::vector<Segment_wrapper_2>& wraps) const {

      Point_2 p, q;
      for (auto& wrap : wraps) {
        const std::size_t group_index = wrap.object_index;
        CGAL_assertion(group_index >= 0 && group_index < lines.size());
        const auto& line = lines[group_index];

        const auto& source = wrap.segment.source();
        const auto& target = wrap.segment.target();

        p = line.projection(source);
        q = line.projection(target);
        wrap.segment = Segment_2(p, q);
      }
    }

    void update_group(
      std::vector<Segment_wrapper_2>& wraps) const {

      Segment_wrapper_2 orth;
      std::vector<Segment_wrapper_2> updated;

      std::size_t count = 0;
      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrap = wraps[i];
        
        // Add a collinear segment.
        wrap.index = count; ++count;
        wrap.is_used = false;
        updated.push_back(wrap);

        // Handle last segment.
        std::size_t j = std::size_t(-1);
        if (i < n - 1) j = i + 1;
        if (j == std::size_t(-1)) break;

        // All intermediate segments.
        const auto& wrapi = wraps[i];
        const auto& wrapj = wraps[j];

        const std::size_t groupi = wrapi.object_index;
        const std::size_t groupj = wrapj.object_index;
        if (groupi != groupj) {
          
          const Line_2 line = Line_2(
            wrapj.segment.source(), wrapj.segment.target());
          const auto source = internal::middle_point_2(
            wrapi.segment.source(), wrapi.segment.target());
          const auto target = line.projection(source);
          orth.segment = Segment_2(source, target);

          // Add an orthogonal segment that connects two collinear groups.
          orth.index = count; ++count;
          orth.is_used = false;
          updated.push_back(orth);
        }
      }
      wraps = updated;
    }

    void create_segments_from_groups(
      std::vector<Segment_wrappers_2>& groups,
      std::vector<Segment_wrapper_2>& wraps) const {

      std::size_t count = 0;
      wraps.clear();
      for (auto& group : groups) {
        for (auto& wrap : group) {
          wrap.index = count;
          wrap.is_used = false;
          wraps.push_back(wrap);
          ++count;
        }
      }
    }

    bool connect_contour(
      const FT max_ordinate_2) {

      bool success = false;
      success = clean_segments(m_wraps);
      if (!success) return false;

      if (m_verbose)
        std::cout << "* number of clean segments = " << m_wraps.size() << std::endl;

      /* make_segments_collinear(max_ordinate_2, m_wraps); */
      intersect_segments(m_wraps);
      return success;
    }

    bool clean_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      // Clean.
      remove_zero_length_segments(wraps);
      if (wraps.size() == 0) return false;
      
      // Filter out wrong segments.
      filter_out_wrong_segments(wraps);
      if (wraps.size() == 0) return false;

      return true;
    }

    void filter_out_wrong_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<Segment_wrapper_2> filtered;
      const std::size_t n = wraps.size();
      const std::size_t start = find_initial_index(wraps);

      std::size_t count = 0;
      std::size_t i = start;
      std::vector<Segment_wrapper_2> parallel;
      do {

        const bool success = get_parallel_segments(
          wraps, parallel, i);
        CGAL_assertion(parallel.size() != 0);
        if (!success) return;

        Segment_2 segment;
        const FT sum_length = 
          parallel_segments_to_segment(parallel, segment);
        if (parallel.size() > 1) {
          
          Segment_wrapper_2 wrap;
          wrap.segment = segment;
          wrap.index = count; ++count;
          filtered.push_back(wrap);

        } else if (parallel.size() == 1) {
          
          auto& wrap = parallel[0];
          wrap.index = count; ++count;
          filtered.push_back(wrap);
        }

      } while (i < n && count <= n);
      if (count > n) return;
      wraps = filtered;
    }

    std::size_t find_initial_index(
      const std::vector<Segment_wrapper_2>& wraps) const {
      return 0;
    }

    bool is_parallel_segment(
      const Segment_2& si, const Segment_2& sp) const {

      const FT angle_pi_2 = CGAL::abs(internal::angle_2_degrees(si, sp));
      const bool target_cond = ( angle_pi_2 <= m_angle_threshold_2 );
      return target_cond;
    }

    bool get_parallel_segments(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<Segment_wrapper_2>& parallel,
      std::size_t& seed) const {
        
      parallel.clear();
      const std::size_t n = wraps.size();
      
      std::size_t i = seed;
      bool next_is_parallel = false;
      std::size_t count = 0;
      do {

        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;
        parallel.push_back(wraps[i]);

        if (ip == std::size_t(-1)) {
          seed = n; return true;
        }

        const auto& si = wraps[i].segment;
        const auto& sp = wraps[ip].segment;
        next_is_parallel = is_parallel_segment(si, sp);
        i = ip;

        ++count;
      } while (next_is_parallel && count < n);
      if (count >= n) return false;
      seed = i;
      return true;
    }

    FT parallel_segments_to_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      Segment_2& result) const {

      Segment_2 ref_segment = find_weighted_segment(wraps);
      const Line_2 line = Line_2(
        ref_segment.source(), ref_segment.target());
      
      FT sum_length = FT(0);
      std::vector<Point_2> points;
      for (const auto& wrap : wraps) {
        
        const Point_2 source = line.projection(wrap.segment.source());
        const Point_2 target = line.projection(wrap.segment.target());

        points.push_back(source);
        points.push_back(target);

        const Segment_2 segment = Segment_2(source, target);
        sum_length += internal::length_2(segment);
      }
      update_segment(points, ref_segment);
      result = ref_segment;
      return sum_length;
    }

    void update_segment(
      const std::vector<Point_2>& points,
      Segment_2& segment) const {

      FT min_proj_value =  internal::max_value<FT>();
      FT max_proj_value = -internal::max_value<FT>();

      const Vector_2 ref_vector = segment.to_vector();
      const Point_2 ref_point = internal::barycenter_2(points);
      
      Point_2 source, target;
      for (const auto& point : points) {
        const Vector_2 curr_vector(ref_point, point);
        const FT value = CGAL::scalar_product(curr_vector, ref_vector);
        
        if (value < min_proj_value) {
          min_proj_value = value;
          source = point; 
        }
        if (value > max_proj_value) {
          max_proj_value = value;
          target = point; 
        }
      }
      segment = Segment_2(source, target);
    }

    void make_segments_collinear(
      const FT max_ordinate_2,
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<std::size_t> seeds;
      seeds.reserve(wraps.size());
      for (std::size_t i = 0; i < wraps.size(); ++i)
        seeds.push_back(i);

      std::vector<Segment_wrappers_2> groups;
      create_collinear_groups(max_ordinate_2, wraps, seeds, groups);

      const std::size_t before = wraps.size();
      std::vector<Line_2> lines;
      create_lines(groups, lines);
      move_segments_toward_lines(lines, wraps);
      const std::size_t after = wraps.size();
      CGAL_assertion(after >= before);

      if (m_verbose)
        std::cout << 
          "* segments before/after = " << before << "/" << after << std::endl;
    }

    void intersect_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        std::size_t im = std::size_t(-1);
        if (i > 0) im = i - 1;
        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;
        
        auto& si = wraps[i].segment;

        if (im == std::size_t(-1) && ip == std::size_t(-1))
          break;

        if (im == std::size_t(-1) && ip != std::size_t(-1)) {
          const auto& sp = wraps[ip].segment;
          intersect_segment(si, sp);
          continue;
        }
        
        if (im != std::size_t(-1) && ip == std::size_t(-1)) {
          const auto& sm = wraps[im].segment;
          intersect_segment(sm, si);
          continue;
        }

        const auto& sm = wraps[im].segment;
        const auto& sp = wraps[ip].segment;
        intersect_segment(sm, si, sp);
      }
    }

    void intersect_segment(
      Segment_2& si, 
      const Segment_2& sp) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(si.source(), si.target());
      const Line_2 line_2 = Line_2(sp.source(), sp.target());
      const bool success = intersect_2(line_1, line_2, target);

      if (!success) target = si.target();
      si = Segment_2(source, target);
    } 

    void intersect_segment(
      const Segment_2& sm, 
      Segment_2& si) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const bool success = intersect_2(line_1, line_2, source);

      if (!success) source = si.source();
      si = Segment_2(source, target);
    } 

    void intersect_segment(
      const Segment_2& sm, 
      Segment_2& si, 
      const Segment_2& sp) const {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const Line_2 line_3 = Line_2(sp.source(), sp.target());

      const bool success1 = intersect_2(line_1, line_2, source);
      const bool success2 = intersect_2(line_2, line_3, target);

      if (!success1) source = si.source();
      if (!success2) target = si.target();

      si = Segment_2(source, target);
    } 

    bool intersect_2(
      const Line_2& line_1, 
      const Line_2& line_2,
      Point_2& in_point) const {
      
      typename std::result_of<Intersect_2(Line_2, Line_2)>::type result 
      = CGAL::intersection(line_1, line_2);
      if (result) {
        if (const Line_2* line = boost::get<Line_2>(&*result)) 
          return false;
        else {
          const Point_2* point = boost::get<Point_2>(&*result);
          in_point = *point; return true;
        }
      }
      return false;
    }

    template<typename OutputIterator>
    void update_input(
      OutputIterator contour) {

      const std::size_t n = m_wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& wrap = m_wraps[i];
        *(++contour) = wrap.segment.source();
      }
      *(++contour) = m_wraps[n - 1].segment.target();
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H
