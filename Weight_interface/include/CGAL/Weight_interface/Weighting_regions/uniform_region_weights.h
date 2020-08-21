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

#ifndef CGAL_GENERALIZED_UNIFORM_REGION_WEIGHTS_H
#define CGAL_GENERALIZED_UNIFORM_REGION_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the uniform cell in 2D.

    This function always returns 1.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \return the computed area.
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the uniform cell in 2D.

    This function always returns 1.

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

    \sa `CGAL::Generalized_weights::uniform_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_area(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& r) {

    const GeomTraits traits;
    return uniform_area(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the uniform cell in 3D.

    This function always returns 1.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \return the computed area.

    \sa `CGAL::Generalized_weights::uniform_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_area(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the uniform cell in 3D.

    This function always returns 1.

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

    \sa `CGAL::Generalized_weights::uniform_area()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_area(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& r) {

    const GeomTraits traits;
    return uniform_area(p, q, r, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief computes area of the uniform cell.

    This function always returns 1.

    \return the computed area.
  */
  double uniform_area() {
    return 1.0;
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_UNIFORM_REGION_WEIGHTS_H
