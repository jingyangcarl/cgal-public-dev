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
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_WEIGHT_INTERFACE_INVERSE_DISTANCE_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_INVERSE_DISTANCE_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace inverse_distance_ns {

    template<typename FT>
    const FT weight(const FT d) {

      FT w = FT(0);
      CGAL_assertion(d != FT(0));
      if (d != FT(0))
        w = FT(1) / d;
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight in 2D or 3D.

    The weight is computed as
    \f$w = \frac{1}{d}\f$
    with notations shown in the figure below.

    This weight is a special case of the `shepard_weight()`.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \cgalFigureBegin{inverse_distance_weight, inverse_distance.svg}
      Notation used for the inverse distance weight.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`

    \param p
    the first point

    \param q
    the second point

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type

    \cgalModels `analytic_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point& q,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight in 2D or 3D.

    This function calls the function `inverse_distance_weight()`.

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`

    \param p
    the first point

    \param q
    the second point

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type

    \sa `inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point& q,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_2(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    const GeomTraits traits;
    return inverse_distance_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return inverse_distance_weight(stub, p, stub, q, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    CGAL::Point_2<GeomTraits> stub;
    return inverse_distance_weight(stub, p, stub, q);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_3(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    const GeomTraits traits;
    return inverse_distance_weight(t, r, p, q, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return inverse_distance_weight(stub, p, stub, q, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    CGAL::Point_3<GeomTraits> stub;
    return inverse_distance_weight(stub, p, stub, q);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_INVERSE_DISTANCE_WEIGHTS_H
