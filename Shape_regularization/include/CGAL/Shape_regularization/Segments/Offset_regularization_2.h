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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_wrapper_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>
#include <CGAL/Shape_regularization/internal/Offset_conditions_2.h>

// TODO:
// * Clean it up.

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief An offset-based regularization type for 2D segments that preserves 
    collinearity relationships.

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
    using Direction_2 = typename Traits::Direction_2;

    using Segment_wrapper_2 = typename internal::Segment_wrapper_2<Traits>;
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

      \tparam NamedParameters
      a sequence of \ref pmp_namedparameters "Named Parameters".

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param np
      optional sequence of \ref pmp_namedparameters "Named Parameters" 
      among the ones listed below

      \param max_offset
      max offset bound in meters, the default is 0.1 meters 

      \param segment_map
      an instance of `SegmentMap` that maps an item from input range to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 1`
      \pre `max_offset >= 0`
    */
    template<typename NamedParameters>
    Offset_regularization_2(
      InputRange& input_range,
      const NamedParameters np,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_num_modified_segments(0) { 

      CGAL_precondition(input_range.size() > 1);
      const FT max_offset = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_offset), FT(1) / FT(10));
      CGAL_precondition(max_offset >= FT(0));

      m_max_offset = max_offset;
      if (m_max_offset < FT(0)) {
        std::cout << "WARNING: The max offset bound has to be within [0, +inf)! ";
        std::cout << " Setting to the default value: 1/10 meters." << std::endl;
        m_max_offset = FT(1) / FT(10);
      }
      clear();
    }
    
    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief calculates the target value between 2 segments, which are
      direct neighbors to each other. The target value is the offset distance.

      \param i
      index of the first segment

      \param j
      index of the second segment

      \pre `i >= 0 && i < input_range.size()`
      \pre `j >= 0 && j < input_range.size()`
    */
    FT target(
      const std::size_t i, 
      const std::size_t j) {

      if (m_wraps.size() == 0) 
        return FT(0);

      CGAL_assertion(m_wraps.size() == m_input_range.size());
      CGAL_assertion(i >= 0 && i < m_wraps.size());
      CGAL_assertion(j >= 0 && j < m_wraps.size());

      const auto& wrapi = m_wraps[i];
      CGAL_assertion(wrapi.is_used);
      const auto& wrapj = m_wraps[j];
      CGAL_assertion(wrapj.is_used);

      const FT tar_val = 
        wrapi.ref_coords.y() - wrapj.ref_coords.y();
      if (CGAL::abs(tar_val) < bound(i) + bound(j))
        m_targets[std::make_pair(i, j)] = tar_val;
      return tar_val;
    }

    /*!
      \brief returns `max_offset`.
    */
    FT bound(const std::size_t) const {
      return m_max_offset;
    }

    /*!
      \brief applies new positions computed by the QP solver 
      to the initial segments.

      \param solution
      a vector with updated segment positions.

      \pre `solution.size() > 0`
    */
    void update(const std::vector<FT>& solution) {
      CGAL_precondition(solution.size() > 0);

      Targets_map targets;
      std::map<FT, Indices> collinear_groups;
      std::map<std::size_t, Segment_wrapper_2> segments;

      CGAL_precondition(m_targets.size() > 0);
      for (const auto& group : m_groups) {
        if (group.size() < 2) continue; 

        targets.clear(); segments.clear();
        build_grouping_data(group, segments, targets);

        collinear_groups.clear();
        if (segments.size() > 0) {
          const std::size_t n = m_input_range.size();

          m_grouping.make_groups(
            m_max_offset, n, segments, solution, 
            collinear_groups, targets);
          translate_collinear_segments(collinear_groups);
        }
      }
    }
    /// @}

    /// \name Miscellaneous
    /// @{ 

    /*!
      \brief inserts a group of segments from `input_range`.

      \tparam IndexRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.
      The value type is `std::size_t`.

      \param index_range
      a const range of segment indices

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
      \brief inserts all input segments from `input_range` as one unique group.

      For more details, 
      see `CGAL::Shape_regularization::Offset_regularization_2::add_group()`.
    */
    void create_unique_group() {
      
      Indices group(m_input_range.size());
      std::iota(group.begin(), group.end(), 0);
      add_group(group);
    }

    /*!
      \brief returns the number of modifed segments.
    */
    std::size_t number_of_modified_segments() const {
      return m_num_modified_segments;
    }

    /// @}

    /// \name Internal data management
    /// @{ 

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_wraps.clear();
      m_wraps.resize(m_input_range.size());
    }

    /// @}

  private:
    Input_range& m_input_range;
    const Segment_map m_segment_map;
    FT m_max_offset;

    std::vector<Segment_wrapper_2> m_wraps;

    std::map<Size_pair, FT> m_targets;
    Grouping m_grouping;
    std::vector<Indices> m_groups;
    std::size_t m_num_modified_segments;

    void update_segment_data(const Indices& group) {
      if (group.size() < 2) return;

      Point_2 frame_origin;
      for (std::size_t i = 0; i < group.size(); ++i) {
        const std::size_t seg_index = group[i];
        CGAL_assertion(
          seg_index >= 0 && seg_index < m_wraps.size());
        
        auto& wrap = m_wraps[seg_index];
        const auto& segment = 
          get(m_segment_map, *(m_input_range.begin() + seg_index));
        wrap.set_all(seg_index, segment);

        if (i == 0)
          frame_origin = wrap.barycenter;
        wrap.set_ref_coords(frame_origin);
      } 
    }

    void build_grouping_data(
      const Indices& group,
      std::map<std::size_t, Segment_wrapper_2>& segments,
      Targets_map& targets) {
      
      for (const std::size_t seg_index : group) {
        CGAL_assertion(
          seg_index >= 0 && seg_index < m_wraps.size());
        const auto& seg_data = m_wraps[seg_index];

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
        CGAL_assertion(
          longest >= 0 && longest < m_wraps.size());
        const auto& longest_data = m_wraps[longest];

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
            CGAL_assertion(
              seg_index >= 0 && seg_index < m_wraps.size());
            const auto& seg_data = m_wraps[seg_index];

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
        const FT seg_length = m_wraps[seg_index].length;
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
      auto& seg_data = m_wraps[seg_index];

      const auto& direction = seg_data.direction;
      const Vector_2 final_normal = Vector_2(
        -direction.dy(), direction.dx());

      const auto& segment = get(m_segment_map, 
        *(m_input_range.begin() + seg_index));

      const auto& source = segment.source();
      const auto& target = segment.target();

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
      const Direction_2& direction) {
      
      FT difference = new_difference;
      auto& seg_data = m_wraps[seg_index];

      seg_data.direction = direction;
      if (seg_data.direction.dy() < FT(0) || 
      (seg_data.direction.dy() == FT(0) && seg_data.direction.dx() < FT(0))) 
        seg_data.direction = -seg_data.direction;

      Vector_2 final_normal = Vector_2(
        -seg_data.direction.dy(), seg_data.direction.dx());

      const auto& segment = get(m_segment_map, 
        *(m_input_range.begin() + seg_index));

      const auto& source = segment.source();
      const auto& target = segment.target();

      FT x1, x2, y1, y2;
      if (
        CGAL::abs(seg_data.direction.dx()) > 
        CGAL::abs(seg_data.direction.dy())) {

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
