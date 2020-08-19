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

#ifndef CGAL_GENERALIZED_COTANGENT_WEIGHTS_H
#define CGAL_GENERALIZED_COTANGENT_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

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
    must be `FieldNumberType`.

    \param cot
    the cotangent value

    \return the computed half weight.

    \sa `CGAL::Generalized_weights::cotangent_weight_2()`
    \sa `CGAL::Generalized_weights::cotangent_weight_3()`
  */
  template<typename FT>
  const FT half_cotangent_weight(
    const FT cot) {

    return cotangent_ns::half_weight(cot);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the cotangent weight for 2D points.

    The weight is computed as
    \f$w = 2 (\cot\beta + \cot\gamma)\f$
    with notations shown in the figure below.

    This weight is equal to the `CGAL::Generalized_weights::discrete_harmonic_weight_2()`.

    This weight is a special case of the `CGAL::Generalized_weights::three_point_family_weight_2()`.

    \cgalFigureBegin{cotangent_weight_2, cotangent.svg}
      Notation used for the cotangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

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
  decltype(auto) cotangent_weight_2(
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

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the cotangent weight for 2D points.

    The weight is computed as
    \f$w = 2 (\cot\beta + \cot\gamma)\f$
    with notations shown in \cgalFigureRef{cotangent_weight_2}.

    This weight is equal to the `CGAL::Generalized_weights::discrete_harmonic_weight_2()`.

    This weight is a special case of the `CGAL::Generalized_weights::three_point_family_weight_2()`.

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
  decltype(auto) cotangent_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return cotangent_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the cotangent weight for 3D points.

    The weight is computed as
    \f$w = 2 (\cot\beta + \cot\gamma)\f$
    with notations shown in the figure below.

    This weight is equal to the `CGAL::Generalized_weights::discrete_harmonic_weight_3()`.

    This weight is a special case of the `CGAL::Generalized_weights::three_point_family_weight_3()`.

    \cgalFigureBegin{cotangent_weight_3, cotangent.svg}
      Notation used for the cotangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_3`.

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
  decltype(auto) cotangent_weight_3(
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

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the cotangent weight for 3D points.

    The weight is computed as
    \f$w = 2 (\cot\beta + \cot\gamma)\f$
    with notations shown in \cgalFigureRef{cotangent_weight_3}.

    This weight is equal to the `CGAL::Generalized_weights::discrete_harmonic_weight_3()`.

    This weight is a special case of the `CGAL::Generalized_weights::three_point_family_weight_3()`.

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
  decltype(auto) cotangent_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return cotangent_weight_3(q, t, r, p, traits);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_COTANGENT_WEIGHTS_H
