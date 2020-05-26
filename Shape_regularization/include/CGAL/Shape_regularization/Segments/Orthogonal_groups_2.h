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

#ifndef CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H
#define CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H

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
    orthogonal segments.

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
  class Orthogonal_groups_2 {

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
        \cgalParamBegin{max_angle}
          max angle deviation in degrees between two segments, the default is 5 degrees
        \cgalParamEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() > 0`
      \pre `max_angle >= 0 && max_angle <= 90`
    */
    template<typename NamedParameters>
    Orthogonal_groups_2(
      const InputRange& input_range,
      const NamedParameters np,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_grouping(
      input_range, np, segment_map) {

      CGAL_precondition(input_range.size() > 0);
      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(5));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));
      m_max_angle = max_angle;
      make_orthogonal_groups();
    }

    /// @}

    // \name Access
    /// @{

    /*!
      \brief returns indices of orthogonal segments organized into groups.

      \tparam OutputIterator
      must be a model of `OutputIterator`

      \param groups
      an instance of OutputIterator,
      whose value type is `std::vector<std::size_t>`

      \return an output iterator
    */
    template<typename OutputIterator>
    OutputIterator groups(OutputIterator groups) const {
      for (const auto& orthogonal_group : m_orthogonal_groups) {
        const auto& group = orthogonal_group;
        *(groups++) = group;
      }
      return groups;
    }
    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    const Parallel_groups_2 m_grouping;

    FT m_max_angle;
    std::vector<Indices> m_orthogonal_groups;

    void make_orthogonal_groups() {

      std::vector<Indices> parallel_groups;
      m_grouping.groups(
        std::back_inserter(parallel_groups));

      Indices orthogonal_group;
      std::vector<bool> states(parallel_groups.size(), false);
      for (std::size_t i = 0; i < parallel_groups.size(); ++i) {
        if (states[i]) continue;

        CGAL_assertion(parallel_groups[i].size() > 0);
        const std::size_t si_index = parallel_groups[i][0];
        const auto& si = get(
          m_segment_map, *(m_input_range.begin() + si_index));

        states[i] = true;
        orthogonal_group.clear();
        for (const std::size_t seg_index : parallel_groups[i])
          orthogonal_group.push_back(seg_index);

        traverse_group(i, si, parallel_groups,
        states, orthogonal_group);
        m_orthogonal_groups.push_back(orthogonal_group);
      }
      CGAL_assertion(
        m_orthogonal_groups.size() <= parallel_groups.size());
    }

    void traverse_group(
      const std::size_t i,
      const Segment_2& si,
      const std::vector<Indices>& parallel_groups,
      std::vector<bool>& states,
      Indices& orthogonal_group) const {

      for (std::size_t j = i + 1; j < parallel_groups.size(); ++j) {
        if (states[j]) continue;

        CGAL_assertion(parallel_groups[j].size() > 0);
        const std::size_t sj_index = parallel_groups[j][0];
        const auto& sj = get(
          m_segment_map, *(m_input_range.begin() + sj_index));

        const FT angle_2 =
          CGAL::abs(internal::angle_2(si, sj));
        if (angle_2 <= m_max_angle ||
            angle_2 >= FT(90) - m_max_angle) {

          states[j] = true;
          for (const std::size_t seg_index : parallel_groups[j])
            orthogonal_group.push_back(seg_index);
        }
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ORTHOGONAL_GROUPS_2_H
