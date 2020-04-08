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

  template<
  typename GeomTraits,
  typename DirectionRange,
  typename InputRange,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class User_defined_principal_directions_2 {

  public:
    using Traits = GeomTraits;
    using Direction_range = DirectionRange;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;

    using FT_pair = std::pair<FT, FT>;
    using Base = internal::Contour_regularization_base_2<Traits>;
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;

    User_defined_principal_directions_2(
      const DirectionRange& directions,
      const InputRange& input_range,
      const PointMap point_map = PointMap()) :
    m_directions(directions),
    m_input_range(input_range),
    m_point_map(point_map) { 
      CGAL_precondition(directions.size() > 0);
    }

    void estimate(
      const bool is_closed,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (is_closed)
        estimate_closed(bounds, directions, assigned);
      else 
        estimate_open(bounds, directions, assigned);
    }

  private:
    const Direction_range& m_directions;
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Base m_base;

    // Should I add readjust_directions() here as in the multiple class?
    void estimate_closed(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (m_directions.size() == 0)
        return;

      CGAL_assertion(
        m_directions.size() == m_input_range.size());
      std::vector<Segment_wrapper_2> wraps;
      const FT max_value = internal::max_value<FT>();
      m_base.initialize_closed(
        max_value, m_input_range, m_point_map, wraps);

      directions = m_directions;
      bounds.clear(); bounds.reserve(directions.size());
      for (std::size_t i = 0; i < directions.size(); ++i)
        bounds.push_back(std::make_pair(FT(45), FT(45)));
      
      m_base.unify_along_contours_closed(wraps, assigned);
      m_base.correct_directions_closed(wraps, assigned);

      if (m_base.verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : assigned)
          std::cout << direction_index << " ";
        std::cout << std::endl;
      }
    }

    // Should I add readjust_directions() here as in the multiple class?
    void estimate_open(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (m_directions.size() == 0)
        return;

      CGAL_assertion(
        m_directions.size() == m_input_range.size() - 1);
      std::vector<Segment_wrapper_2> wraps;
      const FT max_value = internal::max_value<FT>();
      m_base.initialize_open(
        max_value, m_input_range, m_point_map, wraps);

      directions = m_directions;
      bounds.clear(); bounds.reserve(directions.size());
      for (std::size_t i = 0; i < directions.size(); ++i)
        bounds.push_back(std::make_pair(FT(45), FT(45)));
      
      m_base.unify_along_contours_open(wraps, assigned);
      m_base.correct_directions_open(wraps, assigned);

      if (m_base.verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : assigned)
          std::cout << direction_index << " ";
        std::cout << std::endl;
      }
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_USER_DEFINED_PRINCIPAL_DIRECTIONS_2_H
