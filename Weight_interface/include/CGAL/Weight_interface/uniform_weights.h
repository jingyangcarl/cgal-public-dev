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

#ifndef CGAL_WEIGHT_INTERFACE_UNIFORM_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_UNIFORM_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  #if defined(DOXYGEN_RUNNING)

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the uniform weight in 2D or 3D.

    This function always returns 1.

    The type `GeomTraits::Point` must be either
    `GeomTraits::Point_2` or `GeomTraits::Point_3`.

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2` or `AnalyticWeightTraits_3`

    \cgalModels `analytic_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point&,
    const typename GeomTraits::Point&,
    const GeomTraits&) { }

  #endif // DOXYGEN_RUNNING

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const typename GeomTraits::Point_2&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p) {

    const GeomTraits traits;
    return uniform_weight(q, t, r, p, traits);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_weight(
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const typename GeomTraits::Point_3&,
    const GeomTraits&) {

    using FT = typename GeomTraits::FT;
    return FT(1);
  }

  template<typename GeomTraits>
  const typename GeomTraits::FT uniform_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p) {

    const GeomTraits traits;
    return uniform_weight(q, t, r, p, traits);
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the uniform weight.

    This function always returns 1.
  */
  double uniform_weight() {
    return 1.0;
  }

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_UNIFORM_WEIGHTS_H
