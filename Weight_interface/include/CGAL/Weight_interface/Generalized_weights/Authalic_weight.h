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

#ifndef CGAL_GENERALIZED_AUTHALIC_WEIGHT_H
#define CGAL_GENERALIZED_AUTHALIC_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief Authalic weight.

    The full weight is computed as

    \f$w = 2 \frac{\cot\beta + \cot\gamma}{r^2}\f$

    and the half weight as

    \f$h = 2 \frac{\cot\beta}{r^2}\f$

    with notations shown in the figure below. This weight is equal to the
    `CGAL::Generalized_weights::Wachspress_weight`. This weight is a special
    case of the `CGAL::Generalized_weights::Three_point_family_weight`.

    \cgalFigureBegin{authalic_weight, authalic.svg}
      Notation used for the authalic weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Authalic_weight {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using GT = GeomTraits;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// 2D point type.
    typedef typename GeomTraits::Point_2 Point_2;

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
    Authalic_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief to be added
    */
    const FT cotangent(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return cotangent_2(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT squared_distance(
      const Point_2& p,
      const Point_2& q) const {

      return squared_distance_2(p, q);
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
      \brief to be added
    */
    const FT squared_distance(
      const Point_3& p,
      const Point_3& q) const {

      return squared_distance_3(p, q);
    }

    /*!
      \brief computes the half of the authalic weight.
    */
    const FT operator()(
      const FT cot, const FT d2) const {

      return half_weight(cot, d2);
    }

    /*!
      \brief computes the authalic weight.
    */
    const FT operator()(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      return weight_2(q, t, r, p);
    }

    /*!
      \brief computes the authalic weight.
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

    const FT cotangent_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return internal::cotangent_2(m_traits, p, q, r);
    }

    const FT squared_distance_2(
      const Point_2& p,
      const Point_2& q) const {

      const auto squared_distance_2 =
        m_traits.compute_squared_distance_2_object();
      return squared_distance_2(p, q);
    }

    const FT cotangent_3(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return internal::cotangent_3(m_traits, p, q, r);
    }

    const FT squared_distance_3(
      const Point_3& p,
      const Point_3& q) const {

      const auto squared_distance_3 =
        m_traits.compute_squared_distance_3_object();
      return squared_distance_3(p, q);
    }

    const FT half_weight(
      const FT cot, const FT r2) const {

      FT w = FT(0);
      CGAL_assertion(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = cot * inv;
      }
      return w;
    }

    const FT weight_2(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      const FT cot_gamma = internal::cotangent_2(m_traits, t, r, q);
      const FT cot_beta  = internal::cotangent_2(m_traits, q, r, p);

      const auto squared_distance_2 =
        m_traits.compute_squared_distance_2_object();
      const FT d2 = squared_distance_2(q, r);

      return weight(
        cot_gamma, cot_beta, d2);
    }

    const FT weight_3(
      const Point_3& q,
      const Point_3& t,
      const Point_3& r,
      const Point_3& p) const {

      const FT cot_gamma = internal::cotangent_3(m_traits, t, r, q);
      const FT cot_beta  = internal::cotangent_3(m_traits, q, r, p);

      const auto squared_distance_3 =
        m_traits.compute_squared_distance_3_object();
      const FT d2 = squared_distance_3(q, r);

      return weight(
        cot_gamma, cot_beta, d2);
    }

    const FT weight(
      const FT cot_gamma,
      const FT cot_beta,
      const FT r2) const {

      FT w = FT(0);
      CGAL_assertion(r2 != FT(0));
      if (r2 != FT(0)) {
        const FT inv = FT(2) / r2;
        w = (cot_gamma + cot_beta) * inv;
      }
      return w;
    }
  };

  template<typename Point_2>
  decltype(auto) authalic_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    Authalic_weight<Traits> authalic;
    return authalic(q, t, r, p);
  }

  template<typename Point_3>
  decltype(auto) authalic_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    Authalic_weight<Traits> authalic;
    return authalic(q, t, r, p);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_AUTHALIC_WEIGHT_H
