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

#ifndef CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/Segments/Parallel_groups_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief Organizes segments with a similar orientation into groups of 
    collinear segments.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a model of `ReadablePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Collinear_groups_2 {

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
    using Line_2 = typename Traits::Line_2;
    using Indices = std::vector<std::size_t>;
    using Parallel_groups_2 = Parallel_groups_2<Traits, Input_range, Segment_map>;
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
      an instance of `SegmentMap` that maps an item from `input_range` to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \cgalNamedParamsBegin
        \cgalParamBegin{max_offset} 
          max offset deviation in meters between two segments, the default is 0.2 meters
        \cgalParamEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `max_offset >= 0`
    */
    template<typename NamedParameters>
    Collinear_groups_2(
      const InputRange& input_range,
      const NamedParameters np, 
      const SegmentMap segment_map = SegmentMap()) : 
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(
      input_range, np, segment_map) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_offset = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_offset), FT(1) / FT(5));
      CGAL_precondition(max_offset >= FT(0));
      m_max_offset = max_offset;
      make_collinear_groups();
    }

    /// @}

    // \name Access
    /// @{ 

    /*!
      \brief returns indices of collinear segments organized into groups.

      \tparam OutputIterator 
      must be a model of `OutputIterator`

      \param groups
      an instance of OutputIterator, 
      whose value type is `std::vector<std::size_t>`

      \return an output iterator
    */
    template<typename OutputIterator>
    OutputIterator groups(OutputIterator groups) const {
      for (const auto& collinear_group : m_collinear_groups) {
        const auto& group = collinear_group;
        *(groups++) = group;
      }
      return groups;
    }
    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const Parallel_groups_2 m_grouping;
    
    FT m_max_offset;
    std::vector<Indices> m_collinear_groups;

    void make_collinear_groups() {
      
      std::vector<Indices> parallel_groups;
      m_grouping.groups(
        std::back_inserter(parallel_groups));
      m_collinear_groups.reserve(parallel_groups.size());
      
      Indices collinear_group;
      std::vector<bool> states;

      const FT sq_max_dist = m_max_offset * m_max_offset;
      for (const auto& parallel_group : parallel_groups) {
        CGAL_assertion(parallel_group.size() > 0);

        states.clear();
        states.resize(parallel_group.size(), false);
        handle_parallel_group(
          parallel_group, sq_max_dist, 
          states, collinear_group);
      }
      CGAL_assertion(
        m_collinear_groups.size() >= parallel_groups.size());
    }

    void handle_parallel_group(
      const Indices& parallel_group,
      const FT sq_max_dist,
      std::vector<bool>& states,
      Indices& collinear_group) {

      for (std::size_t i = 0; i < parallel_group.size(); ++i) {
        if (states[i]) continue;
        
        const std::size_t si_index = parallel_group[i];
        const auto& si = get(m_segment_map, 
          *(m_input_range.begin() + si_index));
        
        states[i] = true;
        collinear_group.clear();
        collinear_group.push_back(si_index);

        const Line_2 line = Line_2(si.source(), si.target());
        traverse_group(
          i, line, parallel_group, sq_max_dist, 
          states, collinear_group);
        m_collinear_groups.push_back(collinear_group);
      }
    }

    void traverse_group(
      const std::size_t i,
      const Line_2& line,
      const Indices& parallel_group,
      const FT sq_max_dist,
      std::vector<bool>& states,
      Indices& collinear_group) const {

      for (std::size_t j = i + 1; j < parallel_group.size(); ++j) {
        if (states[j]) continue;
        
        const std::size_t sj_index = parallel_group[j];
        const auto& sj = get(m_segment_map, 
          *(m_input_range.begin() + sj_index));
        
        const auto p = internal::middle_point_2(
          sj.source(), sj.target());
        const auto q = line.projection(p);
        
        const FT sq_dist = CGAL::squared_distance(p, q);
        if (sq_dist <= sq_max_dist) {
          states[j] = true;
          collinear_group.push_back(sj_index);
        }
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_COLLINEAR_GROUPS_2_H
