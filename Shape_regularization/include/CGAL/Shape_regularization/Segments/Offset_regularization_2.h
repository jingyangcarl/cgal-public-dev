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
#include <CGAL/Shape_regularization/Segments/Collinear_groups_2.h>

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
    must be a model of `LvaluePropertyMap` whose key type is the value type of the `InputRange`
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

    using Segment_wrapper_2 = typename internal::Segment_wrapper_2<Traits>;
    using Collinear_groups_2 = Collinear_groups_2<Traits, Input_range, Segment_map>;
    using Indices = std::vector<std::size_t>;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref sr_namedparameters "Named Parameters".

      \param input_range
      an instance of `InputRange` with 2D segments

      \param np
      optional sequence of \ref sr_namedparameters "Named Parameters"
      among the ones listed below

      \param segment_map
      an instance of `SegmentMap` that maps an item from input range to `GeomTraits::Segment_2`,
      if not provided, the default is used

      \cgalNamedParamsBegin
        \cgalParamBegin{max_offset}
          max offset deviation in meters between two segments, the default is 0.5 meters
        \cgalParamEnd
      \cgalNamedParamsEnd

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

      CGAL_precondition(
        input_range.size() > 1);
      const FT max_offset = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_offset), FT(1) / FT(2));
      CGAL_precondition(max_offset >= FT(0));

      m_max_offset = max_offset;
      if (m_max_offset < FT(0)) {
        std::cerr << "WARNING: The max offset bound has to be within [0, +inf)! ";
        std::cerr << " Setting to the default value: 0.5 meters." << std::endl;
        m_max_offset = FT(1) / FT(2);
      }
      clear();
      create_unique_group();
    }

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

      if (m_is_first_call) {
        clear();
        m_is_first_call = false;
      }

      if (index_range.size() < 2) return;
      update_wrappers(index_range);
      ++m_num_groups;
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

      CGAL_precondition(i >= 0 && i < m_input_range.size());
      CGAL_precondition(j >= 0 && j < m_input_range.size());
      CGAL_assertion(m_wraps.size() == m_input_range.size());

      const auto& wrapi = m_wraps[i];
      CGAL_assertion(wrapi.is_used);
      const auto& wrapj = m_wraps[j];
      CGAL_assertion(wrapj.is_used);

      const FT target_value =
        wrapi.ref_coords.y() - wrapj.ref_coords.y();
      return target_value;
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
    void update(
      const std::vector<FT>& solution) {

      CGAL_precondition(solution.size() > 0);
      m_num_modified_segments = 0;
      for (auto& wrap : m_wraps) {
        if (!wrap.is_used) continue;

        const std::size_t seg_index = wrap.index;
        const FT difference = -solution[seg_index];

        const auto& direction = wrap.direction;
        const Vector_2 final_normal = Vector_2(
          -direction.dy(), direction.dx());

        const auto& segment = get(m_segment_map,
          *(m_input_range.begin() + seg_index));
        const auto& source = segment.source();
        const auto& target = segment.target();

        const Point_2 new_source = Point_2(
          source.x() + difference * final_normal.x(),
          source.y() + difference * final_normal.y());
        const Point_2 new_target = Point_2(
          target.x() + difference * final_normal.x(),
          target.y() + difference * final_normal.y());

        const FT bx = (new_source.x() + new_target.x()) / FT(2);
        const FT by = (new_source.y() + new_target.y()) / FT(2);
        wrap.barycenter = Point_2(bx, by);
        wrap.c = -wrap.a * bx - wrap.b * by;

        const Segment_2 modified = Segment_2(new_source, new_target);
        put(m_segment_map,
          *(m_input_range.begin() + seg_index), modified);
        ++m_num_modified_segments;
      }
    }

    /// @}

    /// \name Groups
    /// @{

    /*!
      \brief returns indices of collinear segments organized into groups.

      \tparam OutputIterator
      must be a model of `OutputIterator`

      \param groups
      an instance of `OutputIterator`,
      whose value type is `std::vector<std::size_t>`
    */
    template<typename OutputIterator>
    OutputIterator collinear_groups(OutputIterator groups) const {

      const Collinear_groups_2 grouping(
        m_input_range,
        CGAL::parameters::max_offset(m_max_offset),
        m_segment_map);
      return grouping.groups(groups);
    }

    /// @}

    /// \name Miscellaneous
    /// @{

    /*!
      \brief returns the number of modifed segments.
    */
    std::size_t number_of_modified_segments() const {
      return m_num_modified_segments;
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_num_groups = 0;
      m_is_first_call = true;
      for (auto& wrap : m_wraps)
        wrap.is_used = false;
    }

    /// @}

    // EXTRA METHODS TO TEST THE CLASS!
    /// \cond SKIP_IN_MANUAL
    const std::size_t number_of_groups() const {
      return m_num_groups;
    }
    /// \endcond

  private:
    Input_range& m_input_range;
    const Segment_map m_segment_map;

    FT m_max_offset;
    std::vector<Segment_wrapper_2> m_wraps;

    std::size_t m_num_modified_segments;
    std::size_t m_num_groups;
    bool m_is_first_call;

    void create_unique_group() {

      CGAL_precondition(m_input_range.size() > 1);
      if (m_input_range.size() < 2) return;

      m_wraps.clear();
      m_wraps.resize(m_input_range.size());

      Indices group(m_input_range.size());
      std::iota(group.begin(), group.end(), 0);
      update_segment_data(group);
      m_num_groups = 1;
    }

    template<typename IndexRange>
    void update_segment_data(
      const IndexRange& index_range) {

      Point_2 frame_origin;
      for (std::size_t i = 0; i < index_range.size(); ++i) {
        const auto seg_index = *(index_range.begin() + i);
        CGAL_assertion(
          seg_index >= 0 && seg_index < m_wraps.size());
        auto& wrap = m_wraps[seg_index];

        const auto& segment =
          get(m_segment_map, *(m_input_range.begin() + seg_index));
        wrap.set_qp(seg_index, segment);

        if (i == 0) frame_origin = wrap.barycenter;
        wrap.set_ref_coords(frame_origin);
      }
    }

    template<typename IndexRange>
    void update_wrappers(
      const IndexRange& index_range) {

      Point_2 frame_origin;
      for (std::size_t i = 0; i < index_range.size(); ++i) {
        const auto seg_index = *(index_range.begin() + i);
        CGAL_assertion(
          seg_index >= 0 && seg_index < m_wraps.size());
        auto& wrap = m_wraps[seg_index];

        wrap.is_used = true;

        if (i == 0) frame_origin = wrap.barycenter;
        wrap.set_ref_coords(frame_origin);
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OFFSET_REGULARIZATION_2_H
