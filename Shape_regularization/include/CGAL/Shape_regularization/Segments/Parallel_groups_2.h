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
// Author(s)     : Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

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
    using Indices = std::vector<std::size_t>;
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

      \param max_angle
      max angle deviation between two segments, the default is 5 degrees

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 0`
      \pre `max_angle >= 0`
    */
    template<typename NamedParameters>
    Parallel_groups_2(
      const InputRange& input_range,
      const NamedParameters np, 
      const SegmentMap segment_map = SegmentMap()) : 
    m_input_range(input_range),
    m_segment_map(segment_map) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(5));
      CGAL_precondition(max_angle >= FT(0));
      m_max_angle = std::floor(CGAL::to_double(max_angle));
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
      for (const auto& parallel_group : m_parallel_groups) {
        const auto& group = parallel_group.second;
        *(groups++) = group;
      }
      return groups;
    }
    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    
    double m_max_angle;
    std::map<std::size_t, Indices> m_parallel_groups;

    void make_parallel_groups() {
      
      m_parallel_groups.clear();
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const auto& segment = get(
          m_segment_map, *(m_input_range.begin() + i));

        const auto direction = internal::direction_2(segment).to_vector(); 
        const FT orientation = internal::orientation_2(direction);
        const double fvalue = std::floor(CGAL::to_double(orientation));

        const std::size_t num = static_cast<std::size_t>(
          std::floor(fvalue / m_max_angle));
        const std::size_t key = static_cast<std::size_t>(
          num * m_max_angle);

        m_parallel_groups[key].push_back(i);
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_PARALLEL_GROUPS_2_H
