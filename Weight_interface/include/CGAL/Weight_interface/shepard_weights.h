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

#ifndef CGAL_WEIGHT_INTERFACE_SHEPARD_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_SHEPARD_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace shepard_ns {

    template<typename GeomTraits>
    const typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT d,
      const typename GeomTraits::FT p) {

      using FT = typename GeomTraits::FT;
      FT w = FT(0);
      CGAL_assertion(d != FT(0));
      if (d != FT(0)) {
        FT denom = d;
        if (p != FT(1))
          denom = internal::power(traits, d, p);
        w = FT(1) / denom;
      }
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the Shepard weight in 2D or 3D.

    The weight is computed as
    \f$w = \frac{1}{d^a}\f$
    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 1\f$ this weight is equal to the `inverse_distance_weight()`.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \cgalFigureBegin{shepard_weight, shepard.svg}
      Notation used for the Shepard weight.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param a
    the power parameter

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the Shepard weight in 2D or 3D.

    This function calls the function `shepard_weight()`.

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param a
    the power parameter

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type

    \sa `shepard_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_2(traits, q, r);
    return shepard_ns::weight(traits, d, a);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return shepard_weight(t, r, p, q, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return shepard_weight(stub, p, stub, q, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    CGAL::Point_2<GeomTraits> stub;
    return shepard_weight(stub, p, stub, q, a);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_3(traits, q, r);
    return shepard_ns::weight(traits, d, a);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return shepard_weight(t, r, p, q, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return shepard_weight(stub, p, stub, q, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT shepard_weight(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    CGAL::Point_3<GeomTraits> stub;
    return shepard_weight(stub, p, stub, q, a);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_SHEPARD_WEIGHTS_H
