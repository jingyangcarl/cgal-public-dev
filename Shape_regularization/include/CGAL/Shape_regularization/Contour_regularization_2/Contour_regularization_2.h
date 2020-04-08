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

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Closed_contour_regularization_2.h>
#include <CGAL/Shape_regularization/internal/Open_contour_regularization_2.h>

namespace CGAL {
namespace Shape_regularization {

  struct CLOSED { };
  struct OPEN { };

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

    \tparam ContourDirections
    must be a model of `ContourDirections`.

    \tparam ContourTag
    must be either `CLOSED` or `OPEN`.
  */
  template<
  typename GeomTraits,
  typename ContourDirections,
  typename ContourTag>
  class Contour_regularization_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Contour_directions = ContourDirections;
    using Contour_tag = ContourTag;
    
    using FT = typename Traits::FT;
    using Regularization = typename std::conditional<
      std::is_same<ContourTag, CLOSED>::value,
      internal::Closed_contour_regularization_2<Traits, Contour_directions>,
      internal::Open_contour_regularization_2<Traits, Contour_directions> >::type;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam InputRange
      must be a model of `ConstRange`.

      \tparam PointMap
      must be an `LvaluePropertyMap` whose key type is the value type of the input 
      range and value type is `GeomTraits::Point_2`.

      \tparam NamedParameters
      a sequence of \ref pmp_namedparameters "Named Parameters".

      \param input_range
      a range of points, which form a contour

      \param directions
      estiamted contour directions  

      \param np
      optional sequence of \ref pmp_namedparameters "Named Parameters" 
      among the ones listed below

      \pre `input_range.size() >= 3` for closed contours
      \pre `input_range.size() >= 2` for open contours
    */
    template<
    typename InputRange,
    typename PointMap,
    typename NamedParameters>
    Contour_regularization_2(
      const ContourDirections& directions,
      const InputRange& input_range,
      const PointMap point_map = PointMap(),
      const NamedParameters np = CGAL::parameters::all_default()) { 

      const FT max_offset_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_offset), FT(1) / FT(2));
      m_regularization = std::make_shared<Regularization>(
        directions, max_offset_2);
      m_regularization->initialize(input_range, point_map);
    }

    /// @}

    /// \name Regularization
    /// @{

    /*!
      \brief executes the contour regularization algorithm.

      This method regularizes the contour with respect 
      to the defined principal directions.

      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_2`.

      \param contour
      an `OutputIterator` with contour points
    */
    template<typename OutputIterator>
    void regularize(
      OutputIterator contour) {
      m_regularization->regularize(
        contour);
    }

    /// @}

    /// \name Miscellaneous
    /// @{

    /*!
      \brief returns the number of principal directions in the contour.
    */
    const std::size_t number_of_principal_directions() const {
      return m_regularization->number_of_principal_directions();
    }

    /// @}

  private:
    std::shared_ptr<Regularization> m_regularization;
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
