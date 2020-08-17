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

#ifndef CGAL_GENERALIZED_AUTHALIC_WEIGHTS_H
#define CGAL_GENERALIZED_AUTHALIC_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  // The full weight is computed as
  // \f$w = 2 \frac{\cot\beta + \cot\gamma}{r^2}\f$
  // and the half weight as
  // \f$h = 2 \frac{\cot\beta}{r^2}\f$
  // with notations shown in the figure below. This weight is equal to the
  // `CGAL::Generalized_weights::Wachspress_weight`. This weight is a special
  // case of the `CGAL::Generalized_weights::Three_point_family_weight`.
  // \cgalFigureBegin{authalic_weight, authalic.svg}
  //   Notation used for the authalic weight.
  // \cgalFigureEnd

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

    \tparam FT
    must be `FieldNumberType`.

    \param cot
    the cotangent value

    \param d2
    the squared distance value

    \return the computed half weight.
  */
  template<typename FT>
  const FT half_authalic_weight(
    const FT cot, const FT d2) {

    return authalic_ns::half_weight(cot, d2);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the authalic weight for 2D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) authalic_weight_2(
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

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the authalic weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.
  */
  template<typename Point_2>
  decltype(auto) authalic_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return authalic_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the authalic weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) authalic_weight_3(
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

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the authalic weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.
  */
  template<typename Point_3>
  decltype(auto) authalic_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return authalic_weight_3(q, t, r, p, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_AUTHALIC_WEIGHTS_H
