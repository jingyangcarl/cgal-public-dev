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

#ifndef CGAL_GENERALIZED_TANGENT_WEIGHT_2_H
#define CGAL_GENERALIZED_TANGENT_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief 2D tangent weight.

    The full weight is computed as

    \f$w = 2 \frac{t_m + t}{r}\f$, where \f$t_m = \frac{A_m}{r r_m + D_m}\f$ and
    \f$t = \frac{A}{r r_p + D_p}\f$

    and the half weight as

    \f$h = 2 \frac{t}{r}\f$

    with notations shown in the figure below and dot products

    \f$D_m = (v_j - q) \cdot (v_m - q)\f$ and
    \f$D_p = (v_j - q) \cdot (v_p - q)\f$.

    This weight is equal to the `CGAL::Generalized_weights::Mean_value_weight_2`.
    This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight_2`.

    \cgalFigureBegin{tangent_weight, tangent.svg}
      Notation used for the tangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Tangent_weight_2 {

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
    Tangent_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the half value of 2D tangent weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      return half_weight_2(query, vj, vp);
    }

    /*!
      \brief computes the half value of 2D tangent weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      return half_weight_3(query, vj, vp);
    }

    /*!
      \brief computes 2D tangent weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vm, vj, vp);
    }

    /*!
      \brief computes 2D tangent weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      return weight_3(query, vm, vj, vp);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    const FT half_weight_2(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      const auto sj = vj - query;
      const auto sp = vp - query;

      const FT rj = internal::length_2(m_traits, sj);
      const FT rp = internal::length_2(m_traits, sp);

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT Aj = FT(2) * area_2(vj, vp, query);

      const auto dot_product_2 =
        m_traits.compute_scalar_product_2_object();
      const FT Dj = dot_product_2(sj, sp);

      return half_weight(rj, rp, Aj, Dj);
    }

    const FT half_weight_3(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      const auto sj = vj - query;
      const auto sp = vp - query;

      const FT rj = internal::length_3(m_traits, sj);
      const FT rp = internal::length_3(m_traits, sp);

      const FT Aj = FT(2) *
        internal::area_3(m_traits, vj, vp, query);

      const auto dot_product_3 =
        m_traits.compute_scalar_product_3_object();
      const FT Dj = dot_product_3(sj, sp);

      return half_weight(rj, rp, Aj, Dj);
    }

    const FT half_weight(
      const FT rj, const FT rp,
      const FT Aj, const FT Dj) const {

      const FT Pj = rj * rp + Dj;

      FT w = FT(0);
      CGAL_assertion(Pj != FT(0));
      if (Pj != FT(0)) {
        const FT invj = FT(1) / Pj;
        const FT tj = Aj * invj;
        CGAL_assertion(rj != FT(0));
        if (rj != FT(0)) {
          const FT inv = FT(2) / rj;
          w = tj * inv;
        }
      }
      return w;
    }

    const FT weight_2(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      const auto sm = vm - query;
      const auto sj = vj - query;
      const auto sp = vp - query;

      const FT rm = internal::length_2(m_traits, sm);
      const FT rj = internal::length_2(m_traits, sj);
      const FT rp = internal::length_2(m_traits, sp);

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT Am = FT(2) * area_2(vm, vj, query);
      const FT Aj = FT(2) * area_2(vj, vp, query);

      const auto dot_product_2 =
        m_traits.compute_scalar_product_2_object();
      const FT Dm = dot_product_2(sm, sj);
      const FT Dj = dot_product_2(sj, sp);

      return weight(
        rm, rj, rp, Am, Aj, Dm, Dj);
    }

    const FT weight_3(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      const auto sm = vm - query;
      const auto sj = vj - query;
      const auto sp = vp - query;

      const FT rm = internal::length_3(m_traits, sm);
      const FT rj = internal::length_3(m_traits, sj);
      const FT rp = internal::length_3(m_traits, sp);

      const FT Am = FT(2) * internal::area_3(m_traits, vm, vj, query);
      const FT Aj = FT(2) * internal::area_3(m_traits, vj, vp, query);

      const auto dot_product_3 =
        m_traits.compute_scalar_product_3_object();
      const FT Dm = dot_product_3(sm, sj);
      const FT Dj = dot_product_3(sj, sp);

      return weight(
        rm, rj, rp, Am, Aj, Dm, Dj);
    }

    const FT weight(
      const FT rm, const FT rj, const FT rp,
      const FT Am, const FT Aj,
      const FT Dm, const FT Dj) const {

      const FT Pm = rm * rj + Dm;
      const FT Pj = rj * rp + Dj;

      FT w = FT(0);
      CGAL_assertion(Pm != FT(0) && Pj != FT(0));
      if (Pm != FT(0) && Pj != FT(0)) {
        const FT invm = FT(1) / Pm;
        const FT invj = FT(1) / Pj;
        const FT tm = Am * invm;
        const FT tj = Aj * invj;
        CGAL_assertion(rj != FT(0));
        if (rj != FT(0)) {
          const FT inv = FT(2) / rj;
          w = (tm + tj) * inv;
        }
      }
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_TANGENT_WEIGHT_2_H
