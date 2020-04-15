// Copyright (c) 2020 GeometryFactory Sarl (France).
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

#ifndef CGAL_SHAPE_REGULARIZATION_USER_DEFINED_PRINCIPAL_DIRECTIONS_2_H
#define CGAL_SHAPE_REGULARIZATION_USER_DEFINED_PRINCIPAL_DIRECTIONS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_regularization_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours
    
    \brief Sets multiple user-defined principal directions of the contour.

    This algorithm finds the best-fit edges of the contour with respect to the 
    user-defined principal directions and sets all other necessary data.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange`.

    \tparam PointMap
    must be a `ReadablePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Point_2>`.

    \cgalModels `ContourDirections`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class User_defined_directions_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;

    using FT_pair = std::pair<FT, FT>;
    using Base = internal::Contour_regularization_base_2<Traits>;
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam DirectionRange
      must be a model of `ConstRange`. The value type is `GeomTraits::Direction_2`.

      \param direction_range
      a const range with user-defined principal directions

      \param input_range
      a const range of points, which form a contour

      \param is_closed 
      indicates wether the contour is closed or open

      \param point_map
      an instance of `PointMap`, if not provided, the default is used

      \pre `direction_range.size() > 0`
      \pre `direction_range.size() == input_range.size()` for closed contours
      \pre `direction_range.size() == input_range.size() - 1` for open contours

      \pre `input_range.size() >= 3` for closed contours
      \pre `input_range.size() >= 2` for open contours
    */
    template<typename DirectionRange>
    User_defined_directions_2(
      const DirectionRange& direction_range,
      const InputRange& input_range,
      const bool is_closed,
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_point_map(point_map) { 

      CGAL_precondition(
        direction_range.size() > 0);

      if (is_closed)
        estimate_closed(
          direction_range, m_bounds, m_directions, m_assigned);
      else 
        estimate_open(
          direction_range, m_bounds, m_directions, m_assigned);
    }

    /// @}

    /// \name Directions
    /// @{

    /*!
      \brief orients a given `segment` with the index `query_index` with respect
      to the best user-defined principal direction.

      \param query_index an index of the `segment` in the input contour, in other words, 
      the segment's source point is the point in the contour with the index `query_index`

      \param segment a segment to be oriented

      \pre `query_index >= 0 && query_index < input_range.size()` for closed contours
      \pre `query_index >= 0 && query_index < input_range.size() - 1` for open contours
    */
    void orient(
      const std::size_t query_index,
      Segment_2& segment) const {

      m_base.apply_rotation_to_segment(
        m_bounds, m_directions, m_assigned,
        query_index, segment);
    }

    /// @}

    /// \name Miscellaneous
    /// @{

    /*!
      \brief returns the number of principal directions of the contour.
    */
    const std::size_t number_of_directions() const {
      return m_directions.size();
    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Base m_base;

    std::vector<FT_pair> m_bounds;
    std::vector<Direction_2> m_directions;
    std::vector<std::size_t> m_assigned;

    const bool verbose() const {
      return m_base.verbose();
    }

    template<typename DirectionRange>
    void estimate_closed(
      const DirectionRange& direction_range,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (direction_range.size() == 0)
        return;

      CGAL_precondition(
        direction_range.size() == m_input_range.size());
      std::vector<Segment_wrapper_2> wraps;
      const FT max_value = internal::max_value<FT>();
      m_base.initialize_closed(
        max_value, m_input_range, m_point_map, wraps);

      directions.clear(); directions.reserve(direction_range.size());
      for (const auto& direction : direction_range)
        directions.push_back(direction);

      bounds.clear(); bounds.reserve(directions.size());
      for (std::size_t i = 0; i < directions.size(); ++i)
        bounds.push_back(std::make_pair(FT(45), FT(45)));
      
      m_base.unify_along_contours_closed(wraps, assigned);
      m_base.correct_directions_closed(wraps, assigned);

      // Do we need that here?
      // m_base.readjust_directions(wraps, assigned, directions);

      if (verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : assigned)
          std::cout << direction_index << " ";
        std::cout << std::endl;
      }
    }

    template<typename DirectionRange>
    void estimate_open(
      const DirectionRange& direction_range,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (direction_range.size() == 0)
        return;

      CGAL_precondition(
        direction_range.size() == m_input_range.size() - 1);
      std::vector<Segment_wrapper_2> wraps;
      const FT max_value = internal::max_value<FT>();
      m_base.initialize_open(
        max_value, m_input_range, m_point_map, wraps);

      directions.clear(); directions.reserve(direction_range.size());
      for (const auto& direction : direction_range)
        directions.push_back(direction);

      bounds.clear(); bounds.reserve(directions.size());
      for (std::size_t i = 0; i < directions.size(); ++i)
        bounds.push_back(std::make_pair(FT(45), FT(45)));
      
      m_base.unify_along_contours_open(wraps, assigned);
      m_base.correct_directions_open(wraps, assigned);

      // Do we need that here?
      // m_base.readjust_directions(wraps, assigned, directions);

      if (verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : assigned)
          std::cout << direction_index << " ";
        std::cout << std::endl;
      }
    }
  };

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_USER_DEFINED_PRINCIPAL_DIRECTIONS_2_H
