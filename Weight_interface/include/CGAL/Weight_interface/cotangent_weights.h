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

#ifndef CGAL_WEIGHT_INTERFACE_COTANGENT_WEIGHTS_H
#define CGAL_WEIGHT_INTERFACE_COTANGENT_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Weights {

  /// \cond SKIP_IN_MANUAL
  namespace cotangent_ns {

    template<typename FT>
    const FT half_weight(
      const FT cot) {

      return FT(2) * cot;
    }

    template<typename FT>
    const FT weight(
      const FT cot_beta, const FT cot_gamma) {

      return FT(2) * (cot_beta + cot_gamma);
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefUtils

    \brief computes the half value of the cotangent weight.

    This function constructs the half of the cotangent weight using the precomputed
    cotangent value.

    The returned value is
    \f$2\textbf{cot}\f$.

    \tparam FT
    a model of `FieldNumberType`.

    \param cot
    the cotangent value

    \return the computed half weight.

    \sa `CGAL::Weights::cotangent_weight()`
  */
  template<typename FT>
  const FT half_cotangent_weight(
    const FT cot) {

    return cotangent_ns::half_weight(cot);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the cotangent weight for 2D points.

    The weight is computed as
    \f$w = 2 (\cot\beta + \cot\gamma)\f$
    with notations shown in the figure below.

    This weight is equal to the `CGAL::Weights::discrete_harmonic_weight()`.

    This weight is a special case of the `CGAL::Weights::three_point_family_weight()`.

    \cgalFigureBegin{cotangent_weight, cotangent.svg}
      Notation used for the cotangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_2`.

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
  const typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_2(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_2(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p) {

    const GeomTraits traits;
    return cotangent_weight(q, t, r, p, traits);
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the cotangent weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    \tparam GeomTraits
    a model of `AnalyticWeightTraits_3`.

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

    \sa `CGAL::Weights::cotangent_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_weight(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT cot_beta  = internal::cotangent_3(traits, q, t, r);
    const FT cot_gamma = internal::cotangent_3(traits, r, p, q);
    return cotangent_ns::weight(cot_beta, cot_gamma);
  }

  /// \cond SKIP_IN_MANUAL
  template<typename GeomTraits>
  const typename GeomTraits::FT cotangent_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p) {

    const GeomTraits traits;
    return cotangent_weight(q, t, r, p, traits);
  }
  /// \endcond

} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_COTANGENT_WEIGHTS_H
