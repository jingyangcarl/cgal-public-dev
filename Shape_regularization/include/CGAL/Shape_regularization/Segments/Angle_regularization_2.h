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

#ifndef CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_wrapper_2.h>
#include <CGAL/Shape_regularization/Segments/Parallel_groups_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief An angle-based regularization type for 2D segments that preserves 
    parallelism and orthogonality relationships.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a `ReadablePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.

    \cgalModels `RegularizationType`
  */
  template<
  typename GeomTraits, 
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Angle_regularization_2 {
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
    using Vector_2 = typename Traits::Vector_2;
    using Direction_2 = typename Traits::Direction_2;

    using Segment_wrapper_2 = typename internal::Segment_wrapper_2<Traits>;
    using Parallel_groups_2 = Parallel_groups_2<Traits, Input_range, Segment_map>;
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
      max angle bound in degrees, the default is 25 degrees 

      \param segment_map
      an instance of `SegmentMap` that maps an item from input range to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 1`
      \pre `max_angle >= 0 && max_angle <= 90`
    */
    template<typename NamedParameters>
    Angle_regularization_2(
      InputRange& input_range, 
      const NamedParameters np,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) { 
      
      CGAL_precondition(input_range.size() > 1);
      const FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(25));
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      m_max_angle = max_angle;
      if (m_max_angle < FT(0) || m_max_angle > FT(90)) {
        std::cout << "WARNING: The max angle bound has to be within [0, 90]! ";
        std::cout << " Setting to the default value: 25 degrees." << std::endl;
        m_max_angle = FT(25);
      }
      clear();
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief calculates the target value between 2 segments, which are
      direct neighbors to each other. The traget value is the angle.

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

      const FT diff_ij = wrapi.orientation - wrapj.orientation;
      const int diff_90 = static_cast<int>(
        std::floor(CGAL::to_double(diff_ij / FT(90))));

      const FT to_lower = FT(90) * (static_cast<FT>(diff_90) + FT(0)) - diff_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(diff_90) + FT(1)) - diff_ij;
      
      const FT abs_lower = CGAL::abs(to_lower);
      const FT abs_upper = CGAL::abs(to_upper);
      const FT angle = abs_lower < abs_upper ? to_lower : to_upper;
      const FT target_value = angle;
      return target_value;
    }

    /*!
      \brief returns `max_angle`.
    */
    const FT bound(const std::size_t) const {
      return m_max_angle;
    }

    /*!
      \brief applies new orientations computed by the QP solver 
      to the initial segments.

      \param solution
      a vector with updated segment orientations.

      \pre `solution.size() > 0`
    */
    void update(const std::vector<FT>& solution) {

      CGAL_precondition(solution.size() > 0);
      for (auto& wrap : m_wraps) {
        if (!wrap.is_used) continue;

        wrap.orientation += solution[wrap.index]; 
        const double angle_rad = internal::radians_2(wrap.orientation);

        const FT x = static_cast<FT>(std::cos(angle_rad));
        const FT y = static_cast<FT>(std::sin(angle_rad));
        Vector_2 v = Vector_2(x, y);
        internal::normalize(v);

        wrap.direction = internal::direction_2(v);
        const Direction_2 orth = Direction_2(
          -wrap.direction.dy(), wrap.direction.dx());
        wrap.a = orth.dx();
        wrap.b = orth.dy();
        wrap.c = -wrap.a * wrap.barycenter.x() - wrap.b * wrap.barycenter.y();
        
        put(m_segment_map, 
          *(m_input_range.begin() + wrap.index), wrap.orient());
        ++m_num_modified_segments;
      }
    }

    /// @}

    /// \name Miscellaneous
    /// @{ 

    /*!
      \brief returns indices of parallel segments organized into groups.

      \param groups
      an instance of OutputIterator
    */
    template<typename OutputIterator>
    OutputIterator parallel_groups(OutputIterator groups) const {

      const Parallel_groups_2 grouping(
        m_input_range, 
        CGAL::parameters::max_angle(m_max_angle), 
        m_segment_map);
      return grouping.parallel_groups(groups);
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
      
      CGAL_precondition(index_range.size() > 1);
      if (index_range.size() < 2) return;
      
      Indices group;
      group.reserve(index_range.size());
      for (const auto seg_index : index_range)
        group.push_back(seg_index);
      update_segment_data(group);
      ++m_num_groups;
    }

    /*!
      \brief inserts all input segments from `input_range` as one unique group.

      For more details, 
      see `CGAL::Shape_regularization::Angle_regularization_2::add_group()`.
    */
    void create_unique_group() {
      
      CGAL_precondition(m_input_range.size() > 1);
      if (m_input_range.size() < 2) return;

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
      m_num_modified_segments = 0;
      m_num_groups = 0;
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
    
    FT m_max_angle;
    std::vector<Segment_wrapper_2> m_wraps;
    
    std::size_t m_num_modified_segments;
    std::size_t m_num_groups;

    void update_segment_data(
      const Indices& group) {
      if (group.size() < 2) return;

      Segment_wrapper_2 wrap;
      for (const std::size_t seg_index : group) {
        CGAL_assertion(
          seg_index >= 0 && seg_index < m_wraps.size());
        auto& wrap = m_wraps[seg_index];
        const auto& segment = get(m_segment_map, 
          *(m_input_range.begin() + seg_index));
        wrap.set_qp(seg_index, segment);
      }
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H
