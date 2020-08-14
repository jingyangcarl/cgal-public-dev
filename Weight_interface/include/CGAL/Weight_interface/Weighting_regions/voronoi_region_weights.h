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

#ifndef CGAL_GENERALIZED_VORONOI_REGION_WEIGHTS_H
#define CGAL_GENERALIZED_VORONOI_REGION_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  // This weight is the area of the shaded region in the figure below. The region
  // is formed by two midpoints of the edges incident to `q` and the circumcenter of
  // the triangle `[p, q, r]`.
  // \cgalFigureBegin{voronoi_region_weight, voronoi_cell.svg}
  //   Notation used for the Voronoi region weight.
  // \cgalFigureEnd

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 2D triangle [p, q, r].

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

    \return the computed area.
  */
  template<typename GeomTraits>
  decltype(auto) voronoi_area_2(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    const auto circumcenter_2 =
      traits.construct_circumcenter_2_object();
    const Point_2 center =
      circumcenter_2(p, q, r);
    const Point_2 m1 =
      internal::barycenter_2(traits, q, r);
    const Point_2 m2 =
      internal::barycenter_2(traits, q, p);

    const FT A1 = internal::positive_area_2(traits, q, m1, center);
    const FT A2 = internal::positive_area_2(traits, q, center, m2);
    return A1 + A2;
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 2D triangle [p, q, r].

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.
  */
  template<typename Point_2>
  decltype(auto) voronoi_area_2(
    const Point_2& p,
    const Point_2& q,
    const Point_2& r) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return voronoi_area_2(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 3D triangle [p, q, r].

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

    \return the computed area.
  */
  template<typename GeomTraits>
  decltype(auto) voronoi_area_3(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;

    const auto circumcenter_3 =
      traits.construct_circumcenter_3_object();
    const Point_3 center =
      circumcenter_3(p, q, r);
    const Point_3 m1 =
      internal::barycenter_3(traits, q, r);
    const Point_3 m2 =
      internal::barycenter_3(traits, q, p);

    const FT A1 = internal::positive_area_3(traits, q, m1, center);
    const FT A2 = internal::positive_area_3(traits, q, center, m2);
    return A1 + A2;
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 3D triangle [p, q, r].

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.
  */
  template<typename Point_3>
  decltype(auto) voronoi_area_3(
    const Point_3& p,
    const Point_3& q,
    const Point_3& r) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return voronoi_area_3(p, q, r, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_VORONOI_REGION_WEIGHTS_H
