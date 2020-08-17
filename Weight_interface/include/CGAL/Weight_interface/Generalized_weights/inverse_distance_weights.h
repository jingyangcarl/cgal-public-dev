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

#ifndef CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHTS_H
#define CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  // The full weight is computed as
  // \f$w = \frac{1}{r}\f$
  // with notations shown in the figure below. This weight is a special case of
  // the `CGAL::Generalized_weights::Shepard_weight`.
  // \cgalFigureBegin{inverse_distance_weight, inverse_distance.svg}
  //   Notation used for the inverse distance weight.
  // \cgalFigureEnd

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

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the inverse distance weight for 2D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param q
    a query point

    \param r
    the neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) inverse_distance_weight_2(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2&,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_2(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the inverse distance weight for 2D points.

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
  decltype(auto) inverse_distance_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return inverse_distance_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the inverse distance weight for 2D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) inverse_distance_weight_2(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return inverse_distance_weight_2(p, stub, q, stub, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the inverse distance weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.
  */
  template<typename Point_2>
  decltype(auto) inverse_distance_weight_2(
    const Point_2& p,
    const Point_2& q) {

    Point_2 stub;
    return inverse_distance_weight_2(p, stub, q, stub);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the inverse distance weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \param q
    a query point

    \param r
    the neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) inverse_distance_weight_3(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3&,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT d =
      internal::distance_3(traits, q, r);
    return inverse_distance_ns::weight(d);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the inverse distance weight for 3D points.

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
  decltype(auto) inverse_distance_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return inverse_distance_weight_3(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the inverse distance weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) inverse_distance_weight_3(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return inverse_distance_weight_3(p, stub, q, stub, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the inverse distance weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.
  */
  template<typename Point_3>
  decltype(auto) inverse_distance_weight_3(
    const Point_3& p,
    const Point_3& q) {

    Point_3 stub;
    return inverse_distance_weight_3(p, stub, q, stub);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHTS_H
