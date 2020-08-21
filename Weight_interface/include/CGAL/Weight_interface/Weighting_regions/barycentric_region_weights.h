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

#ifndef CGAL_GENERALIZED_BARYCENTRIC_REGION_WEIGHTS_H
#define CGAL_GENERALIZED_BARYCENTRIC_REGION_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the barycentric cell in 2D.

    This area is the area of the shaded region in the figure below. The region
    is formed by two midpoints of the edges incident to `q` and the barycenter of
    the triangle `[p, q, r]`.

    \cgalFigureBegin{barycentric_area, barycentric_cell.svg}
      Notation used for the barycentric cell.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

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
  const typename GeomTraits::FT barycentric_area(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    const Point_2 center =
      internal::barycenter_2(traits, p, q, r);
    const Point_2 m1 =
      internal::barycenter_2(traits, q, r);
    const Point_2 m2 =
      internal::barycenter_2(traits, q, p);

    const FT A1 = internal::positive_area_2(traits, q, m1, center);
    const FT A2 = internal::positive_area_2(traits, q, center, m2);
    return A1 + A2;
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the barycentric cell in 2D.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.

    \sa `CGAL::Generalized_weights::barycentric_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT barycentric_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return barycentric_area(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the barycentric cell in 3D.

    This is an overload of the 2D weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \param traits
    an instance of `GeomTraits`

    \return the computed area.

    \sa `CGAL::Generalized_weights::barycentric_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT barycentric_area(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    using Point_3 = typename GeomTraits::Point_3;

    const Point_3 center =
      internal::barycenter_3(traits, p, q, r);
    const Point_3 m1 =
      internal::barycenter_3(traits, q, r);
    const Point_3 m2 =
      internal::barycenter_3(traits, q, p);

    const FT A1 = internal::positive_area_3(traits, q, m1, center);
    const FT A2 = internal::positive_area_3(traits, q, center, m2);
    return A1 + A2;
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the barycentric cell in 3D.

    This is an overload of the 2D weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.

    \sa `CGAL::Generalized_weights::barycentric_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT barycentric_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return barycentric_area(p, q, r, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_BARYCENTRIC_REGION_WEIGHTS_H
