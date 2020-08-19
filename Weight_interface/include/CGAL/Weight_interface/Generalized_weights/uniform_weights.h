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

#ifndef CGAL_GENERALIZED_UNIFORM_WEIGHTS_H
#define CGAL_GENERALIZED_UNIFORM_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the uniform weight for 2D points.

    This function always returns 1.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) uniform_weight_2(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the uniform weight for 2D points.

    This function always returns 1.

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
  decltype(auto) uniform_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return uniform_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the uniform weight.

    This function always returns 1.

    \return the computed weight.
  */
  double uniform_weight_2() {
    return 1.0;
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the uniform weight for 3D points.

    This function always returns 1.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) uniform_weight_3(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the uniform weight for 3D points.

    This function always returns 1.

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
  decltype(auto) uniform_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return uniform_weight_3(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    This function always returns 1.

    \brief computes the uniform weight.

    \return the computed weight.
  */
  double uniform_weight_3() {
    return 1.0;
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_UNIFORM_WEIGHTS_H
