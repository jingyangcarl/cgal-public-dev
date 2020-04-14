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

#ifndef CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>
#include <CGAL/Shape_regularization/internal/Offset_conditions_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief An offset-based regularization type on a set of 2D segments that preserves 
    collinearity relationship.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a `ReadablePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.

    \cgalModels `RegularizationType`
  */
  template<
  typename GeomTraits, 
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Offset_regularization_2 {
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
    using Conditions = typename internal::Offset_conditions_2<Traits>;
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

      \param max_offset
      max distance value in meters, the default is 0.1 meters 

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 1`
      \pre `max_offset >= 0`
    */
    Offset_regularization_2 (
      InputRange& input_range,
      const FT max_offset = FT(1) / FT(10),
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_d_max(CGAL::abs(max_offset)),
    m_segment_map(segment_map),
    m_num_modified_segments(0) { 

      CGAL_precondition(input_range.size() > 1);
      CGAL_precondition(max_offset >= FT(0));

      if (max_offset < FT(0)) {
        std::cout << 
          "WARNING: The max offset bound has to be within [0, +inf)! Setting to 0." 
        << std::endl;
        m_d_max = FT(0);
      }
    }
    
    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief calculates the target value between 2 segments, which are
      direct neighbors to each other. The target value is the distance.

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
      \brief returns `max_offset`.
    */
    FT bound(const std::size_t) const {
      return m_d_max;
    }

    /*!
      \brief applies new positions computed by the QP solver 
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

      \tparam IndexRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

      \param index_range
      an instance of `IndexRange`

      \pre `index_range.size() > 1`
    */
    template<typename IndexRange>
  	void add_group(
      const IndexRange& index_range) { 
      
      CGAL_precondition(index_range.size() > 1);
      if (index_range.size() < 2) return;
      
      Indices group;
      group.reserve(index_range.size());
      for (const auto seg_index : index_range)
        group.push_back(seg_index);
      
      m_groups.push_back(group);
      update_segment_data(group);
    }

    /*!
      \brief returns the number of modifed segments.
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

    void update_segment_data(const Indices& group) {
      if (group.size() < 2) return;

      Point_2 frame_origin;
      for (std::size_t i = 0; i < group.size(); ++i) {
        const std::size_t seg_index = group[i];
        CGAL_precondition(m_segments.find(seg_index) == m_segments.end());
        if (m_segments.find(seg_index) != m_segments.end())
          continue;

        const auto& segment = 
          get(m_segment_map, *(m_input_range.begin() + seg_index));
        Segment_data seg_data(segment, seg_index);

        if (i == 0)
          frame_origin = seg_data.barycenter;

        seg_data.ref_coords = internal::transform_coordinates_2(
          seg_data.barycenter, frame_origin, seg_data.orientation);
        m_segments.emplace(
          seg_index, seg_data);
      } 
    }

    void build_grouping_data(
      const Indices& group,
      std::map<std::size_t, Segment_data>& segments,
      Targets_map& targets) {
      
      for (const std::size_t seg_index : group) {
        CGAL_precondition(m_segments.find(seg_index) != m_segments.end());
        const auto& seg_data = m_segments.at(seg_index);

        segments.emplace(seg_index, seg_data);
        std::size_t tar_index = 0;

        for(const auto& target : m_targets) {
          const std::size_t seg_index_tar_i = target.first.first;
          const std::size_t seg_index_tar_j = target.first.second;
          const FT tar_value = target.second;

          if (seg_index_tar_i == seg_index)
            targets[std::make_pair(seg_index_tar_i, seg_index_tar_j)] = 
              std::make_pair(tar_value, tar_index);
          ++tar_index;
        }
      }
    }

    void translate_collinear_segments(
      const std::map<FT, Indices>& collinear_groups) {
      
      for (const auto& collinear_group : collinear_groups) {
        const FT dt = collinear_group.first;
        const auto& group = collinear_group.second;

        const std::size_t longest = find_longest_segment(group);
        CGAL_assertion(m_segments.find(longest) != m_segments.end());
        const auto& longest_data = m_segments.at(longest);

        FT new_difference = dt - longest_data.ref_coords.y();
        set_difference(longest, new_difference);

        const FT la = longest_data.a;
        const FT lb = longest_data.b;
        const FT lc = longest_data.c;
        const auto& ldirection = longest_data.direction;

        // Translate the other segments, so that they rest 
        // upon the line ax + by + c = 0.
        for (const std::size_t seg_index : group) {
          if (seg_index != longest) {
            CGAL_precondition(m_segments.find(seg_index) != m_segments.end());
            const auto& seg_data = m_segments.at(seg_index);

            new_difference = dt - seg_data.ref_coords.y();
            set_difference(seg_index, new_difference, la, lb, lc, ldirection);
          }
        }
      }
    }

    std::size_t find_longest_segment(
      const Indices& group) const {
      
      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);
      for (const std::size_t seg_index : group) {
        const FT seg_length = m_segments.at(seg_index).length;
        if (max_length < seg_length) {
          longest = seg_index;
          max_length = seg_length; 
        }
      }
      return longest;
    }

    void set_difference(
      const std::size_t seg_index, 
      const FT new_difference) {

      const FT difference = new_difference;
      auto& seg_data = m_segments.at(seg_index);

      const auto& direction = seg_data.direction;
      const Vector_2 final_normal = Vector_2(-direction.y(), direction.x());

      const auto& source = seg_data.segment.source();
      const auto& target = seg_data.segment.target();

      Point_2 new_source = Point_2(
        source.x() + difference * final_normal.x(), 
        source.y() + difference * final_normal.y());
      Point_2 new_target = Point_2(
        target.x() + difference * final_normal.x(), 
        target.y() + difference * final_normal.y());
      
      const FT bx = (new_source.x() + new_target.x()) / FT(2);
      const FT by = (new_source.y() + new_target.y()) / FT(2);

      m_input_range[seg_index] = Segment_2(new_source, new_target);
      seg_data.c = -seg_data.a * bx - seg_data.b * by;
      ++m_num_modified_segments;
    }

    void set_difference(
      const std::size_t seg_index, 
      const FT new_difference, 
      const FT a, const FT b, const FT c, 
      const Vector_2& direction) {
      
      FT difference = new_difference;
      auto& seg_data = m_segments.at(seg_index);

      seg_data.direction = direction;
      if (seg_data.direction.y() < FT(0) || 
      (seg_data.direction.y() == FT(0) && seg_data.direction.x() < FT(0))) 
        seg_data.direction = -seg_data.direction;

      Vector_2 final_normal = Vector_2(
        -seg_data.direction.y(), seg_data.direction.x());

      const auto& source = seg_data.segment.source();
      const auto& target = seg_data.segment.target();

      FT x1, x2, y1, y2;
      if (CGAL::abs(seg_data.direction.x()) > CGAL::abs(seg_data.direction.y())) {
        x1 = source.x() + difference * final_normal.x();
        x2 = target.x() + difference * final_normal.x(); 
        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;
      } else {    
        y1 = source.y() + difference * final_normal.y();
        y2 = target.y() + difference * final_normal.y();
        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }

      const Point_2 new_source = Point_2(x1, y1);
      const Point_2 new_target = Point_2(x2, y2);
      m_input_range[seg_index] = Segment_2(new_source, new_target);
      ++m_num_modified_segments;
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H
