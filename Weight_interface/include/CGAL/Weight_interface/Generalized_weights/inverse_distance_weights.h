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
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 2D points.

    The weight is computed as
    \f$w = \frac{1}{d}\f$
    with notations shown in the figure below.

    This weight is a special case of the `CGAL::Generalized_weights::shepard_weight()`.

    \cgalFigureBegin{inverse_distance_weight, inverse_distance.svg}
      Notation used for the inverse distance weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param q
    a query point

    \param r
    the neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
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
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& q,
    const CGAL::Point_2<GeomTraits>& t,
    const CGAL::Point_2<GeomTraits>& r,
    const CGAL::Point_2<GeomTraits>& p) {

    const GeomTraits traits;
    return inverse_distance_weight(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 2D points.

    This function calls the function `CGAL::Generalized_weights::inverse_distance_weight()`.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_2& p,
    const typename GeomTraits::Point_2& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_2 stub;
    return inverse_distance_weight(p, stub, q, stub, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 2D points.

    This function calls the function `CGAL::Generalized_weights::inverse_distance_weight()`.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_2`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_2<GeomTraits>& p,
    const CGAL::Point_2<GeomTraits>& q) {

    CGAL::Point_2<GeomTraits> stub;
    return inverse_distance_weight(p, stub, q, stub);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param q
    a query point

    \param r
    the neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
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
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& q,
    const CGAL::Point_3<GeomTraits>& t,
    const CGAL::Point_3<GeomTraits>& r,
    const CGAL::Point_3<GeomTraits>& p) {

    const GeomTraits traits;
    return inverse_distance_weight(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    This function calls the function `CGAL::Generalized_weights::inverse_distance_weight()`.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const typename GeomTraits::Point_3& p,
    const typename GeomTraits::Point_3& q,
    const GeomTraits& traits) {

    typename GeomTraits::Point_3 stub;
    return inverse_distance_weight(p, stub, q, stub, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief computes the inverse distance weight for 3D points.

    This is an overload of the 2D weight for 3D points.

    This function calls the function `CGAL::Generalized_weights::inverse_distance_weight()`.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam GeomTraits
    must be a model of `AnalyticWeightTraits_3`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.

    \sa `CGAL::Generalized_weights::inverse_distance_weight()`
  */
  template<typename GeomTraits>
  const typename GeomTraits::FT inverse_distance_weight(
    const CGAL::Point_3<GeomTraits>& p,
    const CGAL::Point_3<GeomTraits>& q) {

    CGAL::Point_3<GeomTraits> stub;
    return inverse_distance_weight(p, stub, q, stub);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHTS_H
