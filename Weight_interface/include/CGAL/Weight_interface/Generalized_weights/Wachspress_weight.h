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

#ifndef CGAL_GENERALIZED_WACHSPRESS_WEIGHT_H
#define CGAL_GENERALIZED_WACHSPRESS_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  // The full weight is computed as

  // \f$w = \frac{C}{A_m A}\f$

  // with notations shown in the figure below. This weight is equal to the
  // `CGAL::Generalized_weights::Authalic_weight`. This weight is a special
  // case of the `CGAL::Generalized_weights::Three_point_family_weight`.

  // \cgalFigureBegin{wachspress_weight, wachspress.svg}
  //   Notation used for the Wachspress weight.
  // \cgalFigureEnd

  /// \cond SKIP_IN_MANUAL
  namespace internal {

    template<typename FT>
    const FT weight(
      const FT A1, const FT A2, const FT C) {

      FT w = FT(0);
      CGAL_assertion(A1 != FT(0) && A2 != FT(0));
      const FT prod = A1 * A2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = C * inv;
      }
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 2D points.

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
  decltype(auto) wachspress_weight_2(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT C  = internal::area_2(traits, t, r, p);
    return internal::weight(A1, A2, C);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 2D points.

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
  decltype(auto) wachspress_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return wachspress_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 3D points.

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
  decltype(auto) wachspress_weight_3(
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
    return wachspress_weight_2(qf, tf, rf, pf, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 3D points.

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
  decltype(auto) wachspress_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return wachspress_weight_3(q, t, r, p, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WACHSPRESS_WEIGHT_H
