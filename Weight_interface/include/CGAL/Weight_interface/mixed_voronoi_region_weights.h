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

#ifndef CGAL_WEIGHT_INTERFACE_MIXED_VORONOI_REGION_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_MIXED_VORONOI_REGION_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the mixed Voronoi cell in 2D or 3D.

    This area is the area of the shaded region in the figure below. The region
    is formed by two midpoints of the edges incident to `q` and the circumcenter of
    the triangle `[p, q, r]`.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \cgalFigureBegin{mixed_voronoi_area, mixed_voronoi_cell.svg}
      Notation used for the mixed Voronoi cell.
    \cgalFigureEnd

    However, unlike the original `voronoi_area()`, if one of the angles in the triangle
    `[p, q, r]` is obtuse and the circumcenter vertex of the region is outside this triangle,
    this vertex is moved to the mid point of the edge `[r, p]` as shown in the figure below.

    \cgalFigureBegin{mixed_voronoi_area_obtuse, mixed_voronoi_cell_obtuse.svg}
      The case with the obtuse angle.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    this parameter can be omitted if the traits class can be deduced from the point type

    \sa `voronoi_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT mixed_voronoi_area(
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point& q,
    const typename GeomTraits::Point& r,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT mixed_voronoi_area(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    const auto angle_2 =
      traits.angle_2_object();
    const auto a1 = angle_2(p, q, r);
    const auto a2 = angle_2(q, r, p);
    const auto a3 = angle_2(r, p, q);

    Point_2 center;
    if (a1 != CGAL::OBTUSE && a2 != CGAL::OBTUSE && a3 != CGAL::OBTUSE) {
      const auto circumcenter_2 =
        traits.construct_circumcenter_2_object();
      center = circumcenter_2(p, q, r);
    } else {
      center = internal::barycenter_2(traits, r, p);
    }

    const Point_2 m1 =
      internal::barycenter_2(traits, q, r);
    const Point_2 m2 =
      internal::barycenter_2(traits, q, p);

    const FT A1 = internal::positive_area_2(traits, q, m1, center);
    const FT A2 = internal::positive_area_2(traits, q, center, m2);
    return A1 + A2;
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT mixed_voronoi_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return mixed_voronoi_area(p, q, r, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT mixed_voronoi_area(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;

    const auto angle_3 =
      traits.angle_3_object();
    const auto a1 = angle_3(p, q, r);
    const auto a2 = angle_3(q, r, p);
    const auto a3 = angle_3(r, p, q);

    Point_3 center;
    if (a1 != CGAL::OBTUSE && a2 != CGAL::OBTUSE && a3 != CGAL::OBTUSE) {
      const auto circumcenter_3 =
        traits.construct_circumcenter_3_object();
      center = circumcenter_3(p, q, r);
    } else {
      center = internal::barycenter_3(traits, r, p);
    }

    const Point_3 m1 =
      internal::barycenter_3(traits, q, r);
    const Point_3 m2 =
      internal::barycenter_3(traits, q, p);

    const FT A1 = internal::positive_area_3(traits, q, m1, center);
    const FT A2 = internal::positive_area_3(traits, q, center, m2);
    return A1 + A2;
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT mixed_voronoi_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return mixed_voronoi_area(p, q, r, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_MIXED_VORONOI_REGION_WEIGHTS_H
