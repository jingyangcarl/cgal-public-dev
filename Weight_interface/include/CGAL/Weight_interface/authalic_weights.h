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

#ifndef CGAL_WEIGHT_INTERFACE_AUTHALIC_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_AUTHALIC_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace authalic_ns {

    template<typename FT>
    const FT half_weight(
      const FT cot, const FT r2) {

      FT w = FT(0);
      CGAL_assertion(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = cot * inv;
      }
      return w;
    }

    template<typename FT>
    const FT weight(
      const FT cot_gamma, const FT cot_beta, const FT r2) {

      FT w = FT(0);
      CGAL_assertion(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = (cot_gamma + cot_beta) * inv;
      }
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefUtils

    \brief computes the half value of the authalic weight.

    This function constructs the half of the authalic weight using the precomputed
    cotangent and squared distance values.

    The returned value is
    \f$\frac{2\textbf{cot}}{\textbf{d2}}\f$.

    \tparam FT
    a model of `FieldNumberType`.

    \param cot
    the cotangent value

    \param d2
    the squared distance value

    \return the computed half weight.

    \sa `CGAL::Weights::authalic_weight()`
  */
  template<typename FT>
  const FT half_authalic_weight(
    const FT cot, const FT d2) {

    return authalic_ns::half_weight(cot, d2);
  }

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the authalic weight for 2D or 3D points.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    The weight is computed as
    \f$w = 2 \frac{\cot\beta + \cot\gamma}{d^2}\f$
    with notations shown in the figure below.

    This weight is equal to the `CGAL::Weights::wachspress_weight()`.

    This weight is a special case of the `CGAL::Weights::three_point_family_weight()`.

    \cgalFigureBegin{authalic_weight, authalic.svg}
      Notation used for the authalic weight.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`.

    \param q
    a query point

    \param p0
    the first neighbor

    \param p1
    the second neighbor

    \param p2
    the third neighbor

    \param traits
    an instance of `GeomTraits`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point& q,
    const typename GeomTraits::Point& p0,
    const typename GeomTraits::Point& p1,
    const typename GeomTraits::Point& p2,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_gamma = internal::cotangent_2(traits, t, r, q);
    const FT cot_beta  = internal::cotangent_2(traits, q, r, p);

    const auto squared_distance_2 =
      traits.compute_squared_distance_2_object();
    const FT d2 = squared_distance_2(q, r);

    return authalic_ns::weight(
      cot_gamma, cot_beta, d2);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT authalic_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p) {

    const GeomTraits traits;
    return authalic_weight(q, t, r, p, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT authalic_weight(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_gamma = internal::cotangent_3(traits, t, r, q);
    const FT cot_beta  = internal::cotangent_3(traits, q, r, p);

    const auto squared_distance_3 =
      traits.compute_squared_distance_3_object();
    const FT d2 = squared_distance_3(q, r);

    return authalic_ns::weight(
      cot_gamma, cot_beta, d2);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT authalic_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p) {

    const GeomTraits traits;
    return authalic_weight(q, t, r, p, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_AUTHALIC_WEIGHTS_H
