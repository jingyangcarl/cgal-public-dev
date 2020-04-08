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

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Aff_transformation_2.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

// Saver.
#include "../../../../examples/Shape_regularization/include/Saver.h"

// TODO:
// * Can I further simplify this class?
// * Can I use squared distance here instead of distance?
// * Improve find_central_segment().
// * Improve orth segments, which are added during an optimization step. They 
// are too far away from the correct place since they are placed in the middle of
// the segment.
// * Do we actually need make_segments_collienar() in the contour connection?
// * Improve intersection.
// * I should use affine transform to rotate all segments.
// * Can I merge closed and open?
// * Can I put all direction related functions to a separate class?

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class Contour_regularization_base_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    struct Segment_wrapper_2 {

      Segment_2 segment;
      Direction_2 direction;
      std::size_t index = std::size_t(-1);
      std::size_t group = std::size_t(-1);
      bool is_valid_direction = false;
      bool is_used = false;
    };

    using FT_pair = std::pair<FT, FT>;
    using Segments_2 = std::vector<Segment_2>;
    using Segment_wrappers_2 = std::vector<Segment_wrapper_2>;

    Contour_regularization_base_2() :
    m_verbose(true),
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

      CGAL_assertion(input_range.size() >= 3);
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
        wrap.direction = 
          internal::segment_to_direction_2(wrap.segment);
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

      CGAL_assertion(input_range.size() >= 2);
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
        wrap.direction = 
          internal::segment_to_direction_2(wrap.segment);
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
        internal::angle_2_degrees(longest, segment));
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
        internal::rotate(angle_deg, directions[k]);
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
        const FT angle = internal::invar90_angle_2_degrees(di, dj);

        angles[direction_index] += angle;
        counts[direction_index] += FT(1);
      }
    }

  private:
    const bool m_verbose;
    const FT m_angle_threshold_2;
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H
