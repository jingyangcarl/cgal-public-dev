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

#ifndef CGAL_GENERALIZED_COTANGENT_WEIGHT_H
#define CGAL_GENERALIZED_COTANGENT_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief Cotangent weight.

    The full weight is computed as

    \f$w = 2 (\cot\beta + \cot\gamma)\f$

    and the half weight as

    \f$h = 2 \cot\gamma\f$

    with notations shown in the figure below. This weight is equal to the
    `CGAL::Generalized_weights::Discrete_harmonic_weight`. This weight is a special
    case of the `CGAL::Generalized_weights::Three_point_family_weight`.

    \cgalFigureBegin{cotangent_weight, cotangent.svg}
      Notation used for the cotangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Cotangent_weight {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using GT = GeomTraits;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    // /// 2D point type.
    // typedef typename GeomTraits::Point_2 Point_2;

    /// 3D point type.
    typedef typename GeomTraits::Point_3 Point_3;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.
    */
    Cotangent_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief to be added
    */
    template<typename Point_2>
    const FT cotangent(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return cotangent_2(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT cotangent(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return cotangent_3(p, q, r);
    }

    /*!
      \brief computes the half of the cotangent weight.
    */
    const FT operator()(
      const FT cot) const {

      return half_weight(cot);
    }

    /*!
      \brief computes the cotangent weight.
    */
    template<typename Point_2>
    const FT operator()(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      return weight_2(q, t, r, p);
    }

    /*!
      \brief computes the cotangent weight.
    */
    const FT operator()(
      const Point_3& q,
      const Point_3& t,
      const Point_3& r,
      const Point_3& p) const {

      return weight_3(q, t, r, p);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    template<typename Point_2>
    const FT cotangent_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return internal::cotangent_2(m_traits, p, q, r);
    }

    const FT cotangent_3(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return internal::cotangent_3(m_traits, p, q, r);
    }

    const FT half_weight(
      const FT cot) const {

      return FT(2) * cot;
    }

    template<typename Point_2>
    const FT weight_2(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      const FT cot_beta  = internal::cotangent_2(m_traits, q, t, r);
      const FT cot_gamma = internal::cotangent_2(m_traits, r, p, q);
      return weight(cot_beta, cot_gamma);
    }

    const FT weight_3(
      const Point_3& q,
      const Point_3& t,
      const Point_3& r,
      const Point_3& p) const {

      const FT cot_beta  = internal::cotangent_3(m_traits, q, t, r);
      const FT cot_gamma = internal::cotangent_3(m_traits, r, p, q);
      return weight(cot_beta, cot_gamma);
    }

    const FT weight(
      const FT cot_beta,
      const FT cot_gamma) const {

      return FT(2) * (cot_beta + cot_gamma);
    }
  };

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the cotangent weight for 2D points.

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
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    const Cotangent_weight<Traits> cotangent;
    return cotangent(q, t, r, p);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the cotangent weight for 3D points.

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
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    const Cotangent_weight<Traits> cotangent;
    return cotangent(q, t, r, p);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_COTANGENT_WEIGHT_H
