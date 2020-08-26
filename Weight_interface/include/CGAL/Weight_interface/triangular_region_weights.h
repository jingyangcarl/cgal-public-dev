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

#ifndef CGAL_WEIGHT_INTERFACE_TRIANGULAR_REGION_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_TRIANGULAR_REGION_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the triangular cell in 2D or 3D.

    This area is the area of the shaded triangle `[p, q, r]` in the figure below.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \cgalFigureBegin{triangular_area, triangular_cell.svg}
      Notation used for the triangular cell.
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
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT triangular_area(
    const typename GeomTraits::Point& p,
    const typename GeomTraits::Point& q,
    const typename GeomTraits::Point& r,
    const GeomTraits& traits) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT triangular_area(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& r,
    const GeomTraits& traits) {

    return internal::positive_area_2(traits, p, q, r);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT triangular_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return triangular_area(p, q, r, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT triangular_area(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& r,
    const GeomTraits& traits) {

    return internal::positive_area_3(traits, p, q, r);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT triangular_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return triangular_area(p, q, r, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_TRIANGULAR_REGION_WEIGHTS_H
