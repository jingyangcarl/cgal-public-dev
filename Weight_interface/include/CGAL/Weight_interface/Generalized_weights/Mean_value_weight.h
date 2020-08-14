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

#ifndef CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H
#define CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  // The full weight is computed as

  // \f$w = \pm 2 \sqrt{\frac{2 (r_m r_p - D)}{(r r_m + D_m)(r r_p + D_p)}}\f$,

  // with notations shown in the figure below and dot products

  // \f$D_m = (v_j - q) \cdot (v_m - q)\f$,
  // \f$D_p = (v_j - q) \cdot (v_p - q)\f$, and
  // \f$D   = (v_m - q) \cdot (v_p - q)\f$.

  // The \f$\pm\f$ sign is a sign of the weight that depends on the configuration.
  // This weight is equal to the `CGAL::Generalized_weights::Tangent_weight`.
  // This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight`.

  // \cgalFigureBegin{mean_value_weight, mean_value.svg}
  //   Notation used for the mean value weight.
  // \cgalFigureEnd

  /// \cond SKIP_IN_MANUAL
  namespace internal {

    template<typename FT>
    const FT sign_of_weight(
      const FT A1, const FT A2, const FT B) {

      if (A1 > FT(0) && A2 > FT(0) && B <= FT(0)) return  FT(1);
      if (A1 < FT(0) && A2 < FT(0) && B >= FT(0)) return -FT(1);
      if (B  > FT(0)) return  FT(1);
      if (B  < FT(0)) return -FT(1);
      return FT(0);
    }

    template<typename GeomTraits>
    const typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT r1,
      const typename GeomTraits::FT r2,
      const typename GeomTraits::FT r3,
      const typename GeomTraits::FT D1,
      const typename GeomTraits::FT D2,
      const typename GeomTraits::FT D,
      const typename GeomTraits::FT sign) {

      using FT = typename GeomTraits::FT;
      using Get_sqrt = internal::Get_sqrt<GeomTraits>;
      const auto sqrt = Get_sqrt::sqrt_object(traits);

      const FT P1 = r1 * r2 + D1;
      const FT P2 = r2 * r3 + D2;

      FT w = FT(0);
      CGAL_assertion(P1 != FT(0) && P2 != FT(0));
      const FT prod = P1 * P2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = FT(2) * (r1 * r3 - D) * inv;
        CGAL_assertion(w >= FT(0));
        w = sqrt(w);
      }
      w *= FT(2); w *= sign;
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the mean value weight for 2D points.

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
  decltype(auto) mean_value_weight_2(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, t);
    const auto v2 = construct_vector_2(q, r);
    const auto v3 = construct_vector_2(q, p);

    const FT l1 = internal::length_2(traits, v1);
    const FT l2 = internal::length_2(traits, v2);
    const FT l3 = internal::length_2(traits, v3);

    const FT D1 = dot_product_2(v1, v2);
    const FT D2 = dot_product_2(v2, v3);
    const FT D  = dot_product_2(v1, v3);

    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT B  = internal::area_2(traits, p, q, t);
    const FT sign = internal::sign_of_weight(A1, A2, B);

    return internal::weight(
      traits, l1, l2, l3, D1, D2, D, sign);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the mean value weight for 2D points.

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
  decltype(auto) mean_value_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return mean_value_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the mean value weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2` and `AnalyticTraits_3`.

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
  decltype(auto) mean_value_weight_3(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const GeomTraits& traits) {

    using Point_2 = typename GeomTraits::Point_2;
    Point_2 qf, tf, rf, pf;
    internal::flatten(
      traits,
      q,  t,  r,  p,
      qf, tf, rf, pf);
    return mean_value_weight_2(qf, tf, rf, pf, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the mean value weight for 3D points.

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
  decltype(auto) mean_value_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return mean_value_weight_3(q, t, r, p, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H
