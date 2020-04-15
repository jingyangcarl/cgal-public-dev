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

#ifndef CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief Organizes segments with a similar orientation into groups of 
    parallel segments.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a `ReadablePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Parallel_groups_2 {

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
    using Segment_2 = typename Traits::Segment_2;
    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Indices = std::vector<std::size_t>;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{
    /*!
      \brief initializes all internal data structures.

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param max_angle
      max angle deviation between two segments, the default is 5 degrees

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 0`
      \pre `max_angle >= 0`
    */
    Parallel_groups_2 (
      const InputRange& input_range, 
      const FT max_angle = FT(1000000),
      const SegmentMap segment_map = SegmentMap()) : 
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_tolerance(CGAL::abs(max_angle)) {

      CGAL_precondition(input_range.size() > 0);
      CGAL_precondition(max_angle > FT(0));

      build_segment_data();
      make_parallel_groups();
    }

    /// @}

    // \name Access
    /// @{ 

    /*!
      \brief returns indices of parallel segments organized into groups.

      \param groups
      an instance of OutputIterator
    */
    template<typename OutputIterator>
    OutputIterator parallel_groups(OutputIterator groups) {
      
      CGAL_precondition(m_parallel_groups.size() > 0);
      for(const auto& parallel_group : m_parallel_groups) {
        const auto& group = parallel_group.second;
        *(groups++) = group;
      }
      return groups;
    }
    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    std::vector<Segment_data> m_segments;
    const FT m_tolerance;
    std::map<FT, Indices> m_parallel_groups;

    void build_segment_data() {
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const auto& segment = 
          get(m_segment_map, *(m_input_range.begin() + i));
        const Segment_data seg_data(segment, i);
        m_segments.push_back(seg_data);
      }
      CGAL_assertion(m_segments.size() > 0);
    }

    void make_parallel_groups() {
      for (const auto& seg_data : m_segments) {
        const FT angle = static_cast<FT>(floor(
          CGAL::to_double(seg_data.orientation * m_tolerance))) / m_tolerance;
        const std::size_t seg_index = seg_data.index;
        m_parallel_groups[angle].push_back(seg_index);
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
