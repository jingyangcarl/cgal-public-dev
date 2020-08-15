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

#ifndef CGAL_GENERALIZED_WEIGHTS_UTILS_H
#define CGAL_GENERALIZED_WEIGHTS_UTILS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {
namespace utils {

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 2D tangent.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    an instance of `GeomTraits`

    \return the computed tangent.
  */
  template<typename GeomTraits>
  decltype(auto) tangent_2(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    return internal::tangent_2(traits, p, q, r);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 2D tangent.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed tangent.
  */
  template<typename Point_2>
  decltype(auto) tangent_2(
    const Point_2& p,
    const Point_2& q,
    const Point_2& r) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return tangent_2(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 3D tangent.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    an instance of `GeomTraits`

    \return the computed tangent.
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT tangent_3(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    return internal::tangent_3(traits, p, q, r);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 3D tangent.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed tangent.
  */
  template<typename Point_3>
  decltype(auto) tangent_3(
    const Point_3& p,
    const Point_3& q,
    const Point_3& r) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return tangent_3(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 2D cotangent.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    an instance of `GeomTraits`

    \return the computed cotangent.
  */
  template<typename GeomTraits>
  decltype(auto) cotangent_2(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    return internal::cotangent_2(traits, p, q, r);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 2D cotangent.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed cotangent.
  */
  template<typename Point_2>
  decltype(auto) cotangent_2(
    const Point_2& p,
    const Point_2& q,
    const Point_2& r) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return cotangent_2(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 3D cotangent.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    an instance of `GeomTraits`

    \return the computed cotangent.
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_3(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    return internal::cotangent_3(traits, p, q, r);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the 3D cotangent.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed cotangent.
  */
  template<typename Point_3>
  decltype(auto) cotangent_3(
    const Point_3& p,
    const Point_3& q,
    const Point_3& r) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return cotangent_3(p, q, r, traits);
  }

  /// \cond SKIP_IN_MANUAL
  template<typename Point_2>
  decltype(auto) squared_distance_2(
    const Point_2& p,
    const Point_2& q) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    const auto squared_distance_2 =
      traits.compute_squared_distance_2_object();
    return squared_distance_2(p, q);
  }

  template<typename Point_3>
  decltype(auto) squared_distance_3(
    const Point_3& p,
    const Point_3& q) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    const auto squared_distance_3 =
      traits.compute_squared_distance_3_object();
    return squared_distance_3(p, q);
  }

  template<typename Point_2>
  decltype(auto) distance_2(
    const Point_2& p,
    const Point_2& q) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return internal::distance_2(traits, p, q);
  }

  template<typename Point_3>
  decltype(auto) distance_3(
    const Point_3& p,
    const Point_3& q) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return internal::distance_3(traits, p, q);
  }

  template<typename Point_2>
  decltype(auto) area_2(
    const Point_2& p,
    const Point_2& q,
    const Point_2& r) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return internal::area_2(traits, p, q, r);
  }

  template<typename Point_3>
  decltype(auto) area_3(
    const Point_3& p,
    const Point_3& q,
    const Point_3& r) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return internal::positive_area_3(traits, p, q, r);
  }

  template<typename Point_2>
  decltype(auto) scalar_product_2(
    const Point_2& p,
    const Point_2& q,
    const Point_2& r) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    const auto scalar_product_2 =
      traits.compute_scalar_product_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, r);
    const auto v2 = construct_vector_2(q, p);
    return scalar_product_2(v1, v2);
  }

  template<typename Point_3>
  decltype(auto) scalar_product_3(
    const Point_3& p,
    const Point_3& q,
    const Point_3& r) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    const auto scalar_product_3 =
      traits.compute_scalar_product_3_object();
    const auto construct_vector_3 =
      traits.construct_vector_3_object();

    const auto v1 = construct_vector_3(q, r);
    const auto v2 = construct_vector_3(q, p);
    return scalar_product_3(v1, v2);
  }
  /// \endcond

} // namespace utils
} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_UTILS_H
