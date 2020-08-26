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

#ifndef CGAL_WEIGHT_INTERFACE_THREE_POINT_FAMILY_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_THREE_POINT_FAMILY_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace three_point_family_ns {

    template<typename GeomTraits>
    const typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT d1,
      const typename GeomTraits::FT d2,
      const typename GeomTraits::FT d3,
      const typename GeomTraits::FT A1,
      const typename GeomTraits::FT A2,
      const typename GeomTraits::FT B,
      const typename GeomTraits::FT p) {

      using FT = typename GeomTraits::FT;
      FT w = FT(0);
      CGAL_assertion(A1 != FT(0) && A2 != FT(0));
      const FT prod = A1 * A2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        FT r1 = d1;
        FT r2 = d2;
        FT r3 = d3;
        if (p != FT(1)) {
          r1 = internal::power(traits, d1, p);
          r2 = internal::power(traits, d2, p);
          r3 = internal::power(traits, d3, p);
        }
        w = (r3 * A1 - r2 * B + r1 * A2) * inv;
      }
      return w;
    }
  }
  /// \endcond

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the three point family weight in 2D or 3D.

    The weight is computed as
    \f$w = \frac{d_2^a A_1 - d^a B + d_1^a A_2}{A_1 A_2}\f$
    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 0\f$ this weight is equal to the
    `wachspress_weight()` and `authalic_weight()`.

    For \f$a = 1\f$ this weight is equal to the
    `mean_value_weight()` and `tangent_weight()`.

    For \f$a = 2\f$ this weight is equal to the
    `discrete_harmonic_weight()` and `cotangent_weight()`.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \cgalFigureBegin{three_point_family_weight, three_point_family.svg}
      Notation used for the three point family weight.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`.

    \param p0
    the first point

    \param p1
    the second point

    \param p2
    the third point

    \param q
    a query point

    \param a
    the power parameter

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point& p0,
    const typename GeomTraits::Point& p1,
    const typename GeomTraits::Point& p2,
    const typename GeomTraits::Point& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d1 = internal::distance_2(traits, q, t);
    const FT d2 = internal::distance_2(traits, q, r);
    const FT d3 = internal::distance_2(traits, q, p);

    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT B  = internal::area_2(traits, p, q, t);

    return three_point_family_ns::weight(
      traits, d1, d2, d3, A1, A2, B, a);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(t, r, p, q, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using Point_2 = typename GeomTraits::Point_2;
    Point_2 tf, rf, pf, qf;
    internal::flatten(
      traits,
      t,  r,  p,  q,
      tf, rf, pf, qf);
    return three_point_family_weight(tf, rf, pf, qf, a, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(t, r, p, q, a, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_THREE_POINT_FAMILY_WEIGHTS_H
