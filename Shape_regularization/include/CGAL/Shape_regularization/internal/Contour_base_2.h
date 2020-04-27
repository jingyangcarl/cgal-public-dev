// Copyright (c) 2020 GeometryFactory SARL (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot
//

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Segment_wrapper_2.h>

// TODO:
// * Simplify this class if possible.
// * Use squared distance here instead of distance.
// * Improve find_central_segment().
// * Improve orth segments, they are too far away.
// * Improve intersection.
// * Use affine transform to rotate all segments if possible.

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class Contour_base_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    using FT_pair = std::pair<FT, FT>;
    using Segment_wrapper_2 = internal::Segment_wrapper_2<Traits>;
    using Segment_wrappers_2 = std::vector<Segment_wrapper_2>;
    using Polyline = std::vector<Point_3>;

    Contour_base_2() :
    m_verbose(false),
    m_angle_threshold_2(FT(5))
    { }

    const bool verbose() const {
      return m_verbose;
    }

    const FT get_angle_threshold_2() const {
      return m_angle_threshold_2;
    }

    // Optimize this one. It doubles input points.
    template<
    typename Input_range,
    typename Point_map>
    void initialize_closed(
      const FT min_length_2,
      const Input_range& input_range,
      const Point_map point_map,
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_precondition(input_range.size() >= 3);
      const std::size_t n = input_range.size();
      
      wraps.clear();
      wraps.reserve(n);
      
      Segment_wrapper_2 wrap;
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;
        
        const auto& source = get(point_map, *(input_range.begin() + i));
        const auto& target = get(point_map, *(input_range.begin() + ip));

        wrap.index = i;
        wrap.segment = Segment_2(source, target);
        auto v = wrap.segment.to_vector();
        wrap.direction = internal::direction_2(v);
        wrap.is_valid_direction = 
          is_valid_principal_direction(min_length_2, wrap.segment);
        wraps.push_back(wrap);
      }
      CGAL_assertion(wraps.size() == n);
    }

    // Optimize this one. It doubles input points.
    template<
    typename Input_range,
    typename Point_map>
    void initialize_open(
      const FT min_length_2,
      const Input_range& input_range,
      const Point_map point_map,
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_precondition(input_range.size() >= 2);
      const std::size_t n = input_range.size();
      
      wraps.clear();
      wraps.reserve(n);
      
      Segment_wrapper_2 wrap;
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;
        
        const auto& source = get(point_map, *(input_range.begin() + i));
        const auto& target = get(point_map, *(input_range.begin() + ip));

        wrap.index = i;
        wrap.segment = Segment_2(source, target);
        auto v = wrap.segment.to_vector();
        wrap.direction = internal::direction_2(v);
        wrap.is_valid_direction = 
          is_valid_principal_direction(min_length_2, wrap.segment);
        wraps.push_back(wrap);
      }
      CGAL_assertion(wraps.size() == n - 1);
    }

    const bool is_valid_principal_direction(
      const FT min_length_2,
      const Segment_2& segment) const {
      
      CGAL_assertion(min_length_2 >= FT(0));
      const FT threshold = min_length_2 * FT(2);
      const FT squared_threshold = threshold * threshold;
      return segment.squared_length() >= squared_threshold;
    }

    void estimate_initial_directions(
      const FT max_angle_2,
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      std::vector<std::size_t> longest_to_short;
      sort_segments_by_length(wraps, longest_to_short);
      CGAL_assertion(longest_to_short.size() == wraps.size());

      bounds.clear(); directions.clear(); assigned.clear();
      assigned.resize(longest_to_short.size(), std::size_t(-1));

      std::size_t group_index = 0;
      std::size_t query_index = std::size_t(-1);
      do {
        query_index = find_next_longest_segment(
          wraps, longest_to_short);
        if (query_index != std::size_t(-1))
          set_next_longest_direction(
            max_angle_2, wraps, query_index, group_index, 
            bounds, directions, assigned);
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
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& longest_to_short) const {

      std::size_t longest = std::size_t(-1);
      for (std::size_t i = 0; i < longest_to_short.size(); ++i) {
        const std::size_t wrap_index = longest_to_short[i];
        const auto& wrap = wraps[wrap_index];

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
      std::vector<Segment_wrapper_2>& wraps,
      const std::size_t query_index,
      const std::size_t group_index,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(query_index != std::size_t(-1));
      CGAL_assertion(group_index != std::size_t(-1));
      
      // Set current longest direction.
      auto& longest = wraps[query_index];
      assigned[query_index] = group_index;
      longest.is_used = true;

      for (auto& wrap : wraps) {
        if (wrap.index == query_index) // skip longest
          continue;

        // Check if another wrap satisifes the conditions.
        if (is_valid_wrap(wrap)) { 
          if (does_satisify_angle_conditions(
            max_angle_2, longest.segment, wrap.segment)) {
            
            assigned[wrap.index] = group_index;
            wrap.is_used = true;
          }
        }
      }

      // Set internals.
      directions.push_back(longest.direction);
      bounds.push_back(std::make_pair(FT(45), FT(45)));
    }

    // Redo this function using a different method for computing angles.
    bool does_satisify_angle_conditions(
      const FT max_angle_2,
      const Segment_2& longest,
      const Segment_2& segment) const {

      CGAL_precondition(
        max_angle_2 >= FT(0) && max_angle_2 <= FT(90));
      const FT bound_min = max_angle_2;
      const FT bound_max = FT(90) - bound_min;

      const FT angle_2 = CGAL::abs(
        internal::angle_2(longest, segment));
      return (angle_2 <= bound_min) || (angle_2 >= bound_max);
    }

    void set_longest_direction(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction(wraps);

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(wraps.size(), 0);
    }

    Direction_2 compute_longest_direction(
      const std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      CGAL_assertion(n != 0);

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& wrap = wraps[i];
        const FT sq_length = wrap.segment.squared_length();
        if (sq_length > max_length) {
          longest = i; max_length = sq_length;
        }
      }

      CGAL_assertion(longest != std::size_t(-1));
      return wraps[longest].direction;
    }

    void unify_along_contours_closed(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(assigned.size() == wraps.size());
      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrap = wraps[i];
        if (wrap.is_used) continue;
        
        std::size_t im = (i + n - 1) % n;
        std::size_t ip = (i + 1) % n;

        bool stop = false;
        std::size_t max_count = 0;
        do {

          if (wraps[im].is_used) {
            assigned[i] = assigned[im];
            wrap.is_used = true;
            break;
          }

          if (wraps[ip].is_used) {
            assigned[i] = assigned[ip]; 
            wrap.is_used = true;
            break;
          }

          im = (im + n - 1) % n;
          ip = (ip + 1) % n;

          if (im == i || ip == i) 
            stop = true;
          ++max_count;

        } while (!stop && max_count < n);
        if (stop || max_count >= n) {
          std::cerr << 
            "Warning: revert back to the first direction!" << std::endl;
          assigned[i] = 0;
        }
      }
    }

    void correct_directions_closed(
      const std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      const std::size_t n = wraps.size();
      std::vector<std::size_t> clean;
      clean.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const std::size_t dm = assigned[im];
        const std::size_t di = assigned[i];
        const std::size_t dp = assigned[ip];

        if (dm != std::size_t(-1) && dm == dp && di != dm)
          clean.push_back(dm);
        else
          clean.push_back(di);
      }
      assigned = clean;
    }

    void unify_along_contours_open(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrap = wraps[i];
        if (wrap.is_used) continue;
        
        std::size_t im = std::size_t(-1);
        if (i > 0) im = i - 1;
        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;

        bool stop = false;
        std::size_t max_count = 0;
        do {

          if (im != std::size_t(-1) && wraps[im].is_used) {
            assigned[i] = assigned[im];
            wrap.is_used = true;
            break;
          }

          if (ip != std::size_t(-1) && wraps[ip].is_used) {
            assigned[i] = assigned[ip]; 
            wrap.is_used = true;
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
          assigned[i] = 0;
        }
      }
    }

    void correct_directions_open(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<std::size_t>& assigned) const {

      const std::size_t n = wraps.size();
      std::vector<std::size_t> clean;
      clean.reserve(n);
      
      for (std::size_t i = 0; i < n; ++i) {
        
        std::size_t im = std::size_t(-1); // fix this, I skip 0 and last
        if (i > 0) im = i - 1;
        std::size_t ip = std::size_t(-1);
        if (i < n - 1) ip = i + 1;

        std::size_t dm = std::size_t(-1);
        if (im != std::size_t(-1)) dm = assigned[im];
        std::size_t di = std::size_t(-1);
        if (i != std::size_t(-1)) di = assigned[i];
        std::size_t dp = std::size_t(-1);
        if (ip != std::size_t(-1)) dp = assigned[ip];

        if (dm != std::size_t(-1) && dm == dp && di != dm)
          clean.push_back(dm);
        else
          clean.push_back(di);
      }
      assigned = clean;
    }

    void readjust_directions(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& assigned,
      std::vector<Direction_2>& directions) const {

      std::vector<FT> angles, counts;
      create_average_angles(wraps, assigned, directions, 
      angles, counts);

      CGAL_assertion(angles.size() == counts.size());
      CGAL_assertion(angles.size() == directions.size());

      for (std::size_t k = 0; k < angles.size(); ++k) {
        CGAL_assertion(counts[k] != FT(0));
        angles[k] /= counts[k];

        const FT angle_deg = angles[k];
        internal::rotate_direction_2(angle_deg, directions[k]);
        
        // Stable version but requires a segment.
        // internal::rotate_segment_2(angle_deg, FT(0), directions[k]);
      }
    }

    void create_average_angles(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::vector<std::size_t>& assigned,
      const std::vector<Direction_2>& directions,
      std::vector<FT>& angles,
      std::vector<FT>& counts) const {

      CGAL_assertion(directions.size() > 0);

      angles.clear();
      angles.resize(directions.size(), FT(0));

      counts.clear();
      counts.resize(directions.size(), FT(0));

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        if (!wrap.is_valid_direction) continue;

        const std::size_t direction_index = assigned[i];
        CGAL_assertion(direction_index != std::size_t(-1));

        const auto& di = directions[direction_index];
        const auto& dj = wrap.direction;

        // Is it necessary here? Can I use angle_2() instead?
        const FT angle = internal::invar90_angle_2(di, dj);

        angles[direction_index] += angle;
        counts[direction_index] += FT(1);
      }
    }

    void apply_rotation_to_segment(
      const std::vector<FT_pair>& bounds,
      const std::vector<Direction_2>& directions,
      const std::vector<std::size_t>& assigned,
      const std::size_t query_index,
      Segment_2& segment) const {

      CGAL_assertion(assigned.size() > 0);
      CGAL_assertion(bounds.size() == directions.size());
      CGAL_assertion(
        query_index >= 0 && 
        query_index < assigned.size());

      const std::size_t direction_index = assigned[query_index];
      if (direction_index == std::size_t(-1))
        return;
       
      CGAL_assertion(
        direction_index >= 0 && 
        direction_index < directions.size());

      const auto& ref_direction = directions[direction_index];
      const auto& ref_bounds = bounds[direction_index];

      auto v = segment.to_vector();
      const Direction_2 seg_direction = 
        internal::direction_2(v);
      rotate_segment(
        ref_bounds, ref_direction, seg_direction, segment);
    }

    void rotate_contour(
      const std::vector<FT_pair>& bounds,
      const std::vector<Direction_2>& directions,
      const std::vector<std::size_t>& assigned,
      std::vector<Segment_wrapper_2>& wraps) const {

      CGAL_assertion(assigned.size() == wraps.size());
      CGAL_assertion(bounds.size() == directions.size());

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const std::size_t direction_index = assigned[i];
        if (direction_index == std::size_t(-1))
          continue;
       
        CGAL_assertion(
          direction_index >= 0 && 
          direction_index < directions.size());

        auto& wrap = wraps[i];
        const auto& ref_direction = directions[direction_index];
        const auto& ref_bounds = bounds[direction_index];

        const auto& seg_direction = wrap.direction;
        rotate_segment(
          ref_bounds, ref_direction, seg_direction, wrap.segment);
      }
    }

    void rotate_segment(
      const FT_pair& bounds,
      const Direction_2& ref_direction, 
      const Direction_2& seg_direction,
      Segment_2& segment) const {

      // Can I use a segment here?
      const FT angle_deg = internal::compute_angle_2(
        ref_direction, seg_direction);
      const FT converted = CGAL::abs(convert_angle_2(angle_deg));
      if (converted <= bounds.first)
        internal::rotate_segment_2(
          angle_deg, FT(180), segment); // parallel case
      if (converted >= bounds.second)
        internal::rotate_segment_2(
          angle_deg, FT(90), segment); // orthogonal case
    }

    void remove_zero_length_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<Segment_wrapper_2> clean;
      for (const auto& wrap : wraps)
        if (wrap.segment.squared_length() > internal::tolerance<FT>())
          clean.push_back(wrap);
      wraps = clean;
    }

    void optimize_group(
      const FT max_offset_2,
      std::vector<Segment_wrapper_2>& wraps) const {

      std::vector<std::size_t> longest_to_short;
      sort_segments_by_length(wraps, longest_to_short);

      // If max_offset_2 = 0 then each segment is in its own group!
      std::vector<Segment_wrappers_2> groups;
      create_collinear_groups(
        max_offset_2, wraps, longest_to_short, groups);

      // if (m_verbose)
      //   std::cout << 
      //     "* number of collinear groups = " << groups.size() << std::endl;

      std::vector<Line_2> lines;
      create_lines(groups, lines);
      move_segments_toward_lines(lines, wraps);
    }

    void create_collinear_groups(
      const FT max_offset_2,
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
        wrapi.group = group_index;
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
              max_offset_2, source, wrapj.segment)) {

              wrapj.is_used = true;
              wrapj.group = group_index;
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
              max_offset_2, source, wrapj.segment)) {

              wrapj.is_used = true;
              wrapj.group = group_index;
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
      const FT max_offset_2,
      const Point_2& source,
      const Segment_2& segment) const {

      CGAL_assertion(max_offset_2 >= FT(0));
      const Line_2 line = Line_2(
        segment.source(), segment.target());
      const auto target = line.projection(source);
      const Segment_2 proj = Segment_2(source, target);
      
      const FT threshold = max_offset_2;
      const FT squared_threshold = threshold * threshold;
      return proj.squared_length() <= squared_threshold;
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
        const FT sq_distance = wrap.segment.squared_length();
        sum_distance += sq_distance;
        weights.push_back(sq_distance);
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
        const std::size_t group_index = wrap.group;
        CGAL_assertion(group_index >= 0 && group_index < lines.size());
        const auto& line = lines[group_index];

        const auto& source = wrap.segment.source();
        const auto& target = wrap.segment.target();

        p = line.projection(source);
        q = line.projection(target);
        wrap.segment = Segment_2(p, q);
      }
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

    std::pair<bool, bool> is_parallel_segment(
      const Segment_2& sm, 
      const Segment_2& si, 
      const Segment_2& sp) const {

      const FT angle_mi_2 = CGAL::abs(internal::angle_2(sm, si));
      const FT angle_pi_2 = CGAL::abs(internal::angle_2(si, sp));

      const bool source_cond = ( angle_mi_2 <= m_angle_threshold_2 );
      const bool target_cond = ( angle_pi_2 <= m_angle_threshold_2 );

      return std::make_pair(source_cond, target_cond);
    }

    bool is_parallel_segment(
      const Segment_2& si, const Segment_2& sp) const {

      const FT angle_pi_2 = CGAL::abs(internal::angle_2(si, sp));
      const bool target_cond = ( angle_pi_2 <= m_angle_threshold_2 );
      return target_cond;
    }

    void parallel_segments_to_segment(
      const std::vector<Segment_wrapper_2>& wraps,
      Segment_2& result) const {

      Segment_2 ref_segment = find_weighted_segment(wraps);
      const Line_2 line = Line_2(
        ref_segment.source(), ref_segment.target());
      
      std::vector<Point_2> points;
      for (const auto& wrap : wraps) {
        
        const Point_2 source = line.projection(wrap.segment.source());
        const Point_2 target = line.projection(wrap.segment.target());

        points.push_back(source);
        points.push_back(target);
      }
      update_segment(points, ref_segment);
      result = ref_segment;
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

    // Do we need this function at all? It does not seem to work well!
    // At least when I call it from the connect_contour().
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

    void export_polylines(
      const std::vector<Segment_wrapper_2>& wraps,
      const std::string file_path) const {

      std::vector<Segment_2> segments;
      segments.reserve(wraps.size());
      for (const auto& wrap : wraps)
        segments.push_back(wrap.segment);
      export_polylines(segments, file_path);
    }

    void export_polylines(
      const std::vector<Segment_2>& segments,
      const std::string file_path) const {
      
      std::vector<Polyline> polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const auto& s = segments[i].source();
        const auto& t = segments[i].target();
        
        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      export_polylines(polylines, file_path);
    }

    void export_polylines(
      const std::vector<Polyline>& polylines,
      const std::string file_path) const {

      if (polylines.size() == 0)
        return;

      std::stringstream out;
      out.precision(20); 

      for (std::size_t i = 0; i < polylines.size(); ++i) {
        const auto &polyline = polylines[i];

        out << polyline.size() << " ";
        for (std::size_t j = 0; j < polyline.size(); ++j)
          out << polyline[j] << " ";
        out << std::endl;
      }
      save(out, file_path + ".polylines");
    }

    void save(
      const std::stringstream& out,
      const std::string path) const {
      
      std::ofstream file(path.c_str(), std::ios_base::out);
      CGAL::set_ascii_mode(file);
      if (!file) {
        std::cout << 
        "Error: cannot save the file: " << path << std::endl; return;
      }
      
      file << out.str() << std::endl; file.close();
      std::cout << 
        "* segments are saved in " << path << std::endl;
    }

  private:
    const bool m_verbose;
    const FT m_angle_threshold_2;
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_BASE_2_H
