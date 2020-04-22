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

#ifndef CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_regularization_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename ContourDirections>
  class Closed_contour_regularization_2 {

  public:
    using Traits = GeomTraits;
    using Contour_directions = ContourDirections;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Direction_2 = typename Traits::Direction_2;
    using Segment_2 = typename Traits::Segment_2;
    using Line_2 = typename Traits::Line_2;

    using Base = internal::Contour_regularization_base_2<Traits>;

    using FT_pair = std::pair<FT, FT>;
    using Segments_2 = std::vector<Segment_2>;
    
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;
    using Segment_wrappers_2 = typename Base::Segment_wrappers_2;

    Closed_contour_regularization_2(
      const Contour_directions& estimator,
      const FT max_offset_2) :
    m_estimator(estimator),
    m_max_offset_2(max_offset_2) 
    { }

    template<
    typename Input_range,
    typename Point_map>
    void initialize(
      const Input_range& input_range,
      const Point_map point_map) {

      const FT max_value = internal::max_value<FT>();
      m_base.initialize_closed(
        max_value, input_range, point_map, m_wraps);
      CGAL_assertion(m_wraps.size() == input_range.size());
    }

    template<typename OutputIterator>
    void regularize(
      OutputIterator contour) {

      CGAL_assertion(m_wraps.size() >= 3);
      if (m_wraps.size() < 3) return;

      rotate_contour(m_wraps);
      if (verbose())
        m_base.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/rotated");

      bool success = optimize_contour(m_max_offset_2, m_wraps);
      if (!success) return;
      if (verbose())
        m_base.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/optimized");

      success = connect_contour(m_wraps);
      if (!success) return;
      if (verbose())
        m_base.export_polylines(
          m_wraps, "/Users/monet/Documents/gsoc/ggr/logs/connected");

      update_input(m_wraps, contour);
    }

  private:
    const Contour_directions& m_estimator;
    const FT m_max_offset_2;
    const Base m_base;

    std::vector<Segment_wrapper_2> m_wraps;
    
    const bool verbose() const {
      return m_base.verbose();
    }

    void rotate_contour(
      std::vector<Segment_wrapper_2>& wraps) const {

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        auto& wrap = wraps[i];
        m_estimator.orient(i, wrap.segment);
      }
    }

    bool optimize_contour(
      const FT max_offset_2,
      std::vector<Segment_wrapper_2>& wraps) const {

      // Clean.
      m_base.remove_zero_length_segments(wraps);
      CGAL_assertion(wraps.size() >= 4);
      if (wraps.size() < 4) return false; // should be at least a quad

      // Create groups of collinear segments.
      std::vector<Segment_wrappers_2> groups;
      create_consecutive_groups(wraps, groups);
      if (verbose())
        std::cout << 
          "* number of consecutive groups = " << groups.size() << std::endl;

      // Optimize all groups with at least two collinear segments.
      std::size_t count = 0;
      for (auto& group : groups) {
        if (group.size() > 1) {
          m_base.optimize_group(max_offset_2, group);
          update_group(group);
          ++count;
        }
      }
      if (verbose())
        std::cout << 
          "* number of optimized groups = " << count << std::endl;

      // Update segments.
      m_base.create_segments_from_groups(groups, wraps);
      if (wraps.size() < 4) return false;
      return true;
    }

    void create_consecutive_groups(
      std::vector<Segment_wrapper_2>& wraps,
      std::vector<Segment_wrappers_2>& groups) const {

      groups.clear();
      for (auto& wrap : wraps)
        wrap.is_used = false;
      std::vector<Segment_wrapper_2> group;

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        auto& wrapi = wraps[i];
        if (wrapi.is_used) continue;
        
        group.clear(); 
        wrapi.is_used = true;
        group.push_back(wrapi);
        
        const std::size_t ip = (i + 1) % n;
        if (ip != 0) {
          for (std::size_t j = ip; j < n; ++j) {
            auto& wrapj = wraps[j];
            const FT angle_2 = CGAL::abs(
              internal::angle_2_degrees(
                wrapi.segment, wrapj.segment));

            if (angle_2 <= m_base.get_angle_threshold_2()) {
              wrapj.is_used = true;
              group.push_back(wrapj); 
            } else break;
          }
        }
        groups.push_back(group);
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
        const std::size_t j = (i + 1) % n;
        if (j == 0) break;

        // All intermediate segments.
        const auto& wrapi = wraps[i];
        const auto& wrapj = wraps[j];

        const std::size_t groupi = wrapi.group;
        const std::size_t groupj = wrapj.group;
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

    bool connect_contour(
      std::vector<Segment_wrapper_2>& wraps) const {
      
      CGAL_assertion(wraps.size() >= 4);
      if (wraps.size() < 4) 
        return false;

      bool success = false;
      success = clean_segments(wraps);
      if (!success) return false;

      if (verbose())
        std::cout << "* number of clean segments = " << 
        wraps.size() << std::endl;

      // Do we need it here?
      // m_base.make_segments_collinear(max_offset_2, wraps);

      intersect_segments(wraps);
      return success;
    }

    bool clean_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      // Clean.
      m_base.remove_zero_length_segments(wraps);
      if (wraps.size() < 4) return false;
      
      // Filter out wrong segments.
      filter_out_wrong_segments(wraps);
      if (wraps.size() < 4) return false;
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
        m_base.parallel_segments_to_segment(parallel, segment);
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

      } while (i != start && count <= n * 2);
      if (count > n * 2) return;
      wraps = filtered;
    }

    std::size_t find_initial_index(
      const std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        const auto& si = wraps[i].segment;
        const auto& sm = wraps[im].segment;
        const auto& sp = wraps[ip].segment;

        const auto pair = m_base.is_parallel_segment(sm, si, sp);
        const bool previous_is_orthogonal = !(pair.first);
        if (previous_is_orthogonal) return i;
      }
      return 0;
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

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& si = wraps[i].segment;
        const auto& sm = wraps[im].segment;
        const auto& sp = wraps[ip].segment;

        parallel.push_back(wraps[i]);
        const auto pair = m_base.is_parallel_segment(sm, si, sp);
        next_is_parallel = pair.second;
        i = ip;

        ++count;
      } while (next_is_parallel && count < n);
      if (count >= n) return false;
      seed = i;
      return true;
    }

    void intersect_segments(
      std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        auto& si = wraps[i].segment;
        const auto& sm = wraps[im].segment;
        const auto& sp = wraps[ip].segment;
        
        m_base.intersect_segment(sm, si, sp);
      }
    }

    template<typename OutputIterator>
    void update_input(
      const std::vector<Segment_wrapper_2>& wraps,
      OutputIterator contour) const {

      for (std::size_t i = 0; i < wraps.size(); ++i) {
        const auto& wrap = wraps[i];
        *(++contour) = wrap.segment.source();
      }
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CLOSED_CONTOUR_REGULARIZATION_2_H
