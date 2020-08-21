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

#ifndef CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHTS_H
#define CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

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

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the three point family weight for 2D points.

    The weight is computed as
    \f$w = \frac{d_2^a A_1 - d^a B + d_1^a A_2}{A_1 A_2}\f$
    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 0\f$ this weight is equal to the
    `CGAL::Generalized_weights::wachspress_weight()` and
    `CGAL::Generalized_weights::authalic_weight()`.

    For \f$a = 1\f$ this weight is equal to the
    `CGAL::Generalized_weights::mean_value_weight()` and
    `CGAL::Generalized_weights::tangent_weight()`.

    For \f$a = 2\f$ this weight is equal to the
    `CGAL::Generalized_weights::discrete_harmonic_weight()` and
    `CGAL::Generalized_weights::cotangent_weight()`.

    \cgalFigureBegin{three_point_family_weight, three_point_family.svg}
      Notation used for the three point family weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param a
    the power parameter

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
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

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the three point family weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param a
    the power parameter, the default is `mean_value_weight()`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::three_point_family_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(q, t, r, p, a, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the three point family weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2` and `AnalyticWeightTraits_3`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param a
    the power parameter

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::three_point_family_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::FT a,
    const GeomTraits& traits) {

    using Point_2 = typename GeomTraits::Point_2;
    Point_2 qf, tf, rf, pf;
    internal::flatten(
      traits,
      q,  t,  r,  p,
      qf, tf, rf, pf);
    return three_point_family_weight(qf, tf, rf, pf, a, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the three point family weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2` and `AnalyticWeightTraits_3`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param a
    the power parameter, the default is `mean_value_weight()`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::three_point_family_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT three_point_family_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p,
    const typename GeomTraits::FT a =
    typename GeomTraits::FT(1)) {

    const GeomTraits traits;
    return three_point_family_weight(q, t, r, p, a, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHTS_H
