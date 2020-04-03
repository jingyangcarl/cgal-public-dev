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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Andreas Fabri
//

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Shape_regularization/enum.h>
#include <CGAL/Shape_regularization/internal/Closed_contour_regularization_2.h>
#include <CGAL/Shape_regularization/internal/Open_contour_regularization_2.h>

// TODO:
// 1. Use named parameters.
// 2. Simplify the design.
// 3. Should I change position of the GeomTraits?
// 4. Should it return the same number of segments as input? Now it is not.
// 5. I think I should let users set arbitrary directions by filling in the m_longest 
// instead of setting indices of other contour edges.
// 6. I think I should let users stop right after the rotation step if they want.

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef
    
    \brief Contour regularization algorithm.

    This algorithm enables to regularize both open and closed contours.

    \tparam GeomTraits 
    must be a model of `Kernel`.
    
    \tparam InputRange 
    must be a model of `ConstRange`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_2`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Contour_regularization_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;

    using Closed_contour = internal::Closed_contour_regularization_2<
      Traits, Input_range, Point_map>;
    using Open_contour = internal::Open_contour_regularization_2<
      Traits, Input_range, Point_map>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range
      a range of points, which form a contour

      \param point_map
      an instance of `PointMap` that maps an item from `input_range` 
      to `GeomTraits::Point_2`

      \param is_closed_contour
      indicates weather the contour is closed or not, 
      the default is `true`

      \pre `input_range.size() >= 3` for closed contour
      \pre `input_range.size() >= 2` for open contour
    */
    Contour_regularization_2(
      Input_range& input_range,
      Point_map point_map,
      const bool is_closed_contour = true) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_is_closed_contour(is_closed_contour),
    m_closed_contour(input_range, point_map),
    m_open_contour(input_range, point_map)
    { }

    /// @}

    /// \name Functions
    /// @{

    /*!
      \brief sets principal directions of the contour.

      This method sets the user-defined principal contour directions. If the 
      direction index is `std::size_t(-1)`, the corresponding contour edge 
      is not regularized.

      \param directions 
      a set of indices where each index represents a contour edge from the `input_range` 
      that is chosen as one of the principal directions

      \pre `directions.size() == input_range.size()` for closed contour
      \pre `directions.size() == input_range.size() - 1` for open contour
    */
    void set_principal_directions(
      const std::vector<std::size_t>& directions) {

      if (m_is_closed_contour)
        m_closed_contour.set_directions(directions);
      else
        m_open_contour.set_directions(directions);
    }

    /*!
      \brief estimates principal directions of the contour.

      This method estimates the principal contour directions automatically.

      \param direction_type
      indicates which type of principal directions should be estimated, 
      the default is `Direction_type::LONGEST`

      \param min_length_2
      a min length of the contour edge in meters that will be set as a principal direction,
      the default is three meters

      \param max_angle_2 
      a max angle difference in degrees between a segment and its collinear or orthogonal orientation,
      the default is 15 degrees

      \pre `min_length_2 > FT(0)`
      \pre `max_angle_2 >= FT(0) && max_angle_2 <= FT(90)`
    */
    void estimate_principal_directions(
      const Direction_type direction_type = Direction_type::LONGEST,
      const FT min_length_2 = FT(3),
      const FT max_angle_2 = FT(15)) {

      if (m_is_closed_contour)
        m_closed_contour.estimate_principal_directions(
          direction_type, min_length_2, max_angle_2);
      else
        m_open_contour.estimate_principal_directions(
          direction_type, min_length_2, max_angle_2);
    }

    /*!
      \brief executes the contour regularization algorithm.

      This method regularizes the contour with respect to the defined principal directions.
      That is it sets all other not principal segments either orthogonal or collinear 
      to the chosen directions.

      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_2`.

      \param contour
      an `OutputIterator` with contour points

      \param max_ordinate_2
      a max distance in meters between two collinear segments that defines if these segments 
      should be merged or not, the default is half a meter.
    */
    template<typename OutputIterator>
    void regularize(
      OutputIterator contour,
      const FT max_ordinate_2 = FT(1) / FT(2)) {
      
      if (m_is_closed_contour)
        m_closed_contour.regularize(contour, max_ordinate_2);
      else
        m_open_contour.regularize(contour, max_ordinate_2);
    }

    /*!
      \brief returns the number of principal directions in the contour.
    */
    std::size_t number_of_principal_directions() const {

      if (m_is_closed_contour)
        return m_closed_contour.number_of_principal_directions();
      else
        return m_open_contour.number_of_principal_directions();
    }

  private:
    Input_range& m_input_range;
    Point_map m_point_map;

    const bool m_is_closed_contour;

    Closed_contour m_closed_contour;
    Open_contour m_open_contour;
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
