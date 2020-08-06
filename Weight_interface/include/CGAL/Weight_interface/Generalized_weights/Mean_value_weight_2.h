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

#ifndef CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H
#define CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief 2D mean value weight.

    The full weight is computed as

    \f$w = \pm 2 \sqrt{\frac{2 (r_m r_p - D)}{(r r_m + D_m)(r r_p + D_p)}}\f$,

    with notations shown in the figure below and dot products

    \f$D_m = (v_j - q) \cdot (v_m - q)\f$,
    \f$D_p = (v_j - q) \cdot (v_p - q)\f$, and
    \f$D = (v_m - q) \cdot (v_p - q)\f$.

    The \f$\pm\f$ sign is a sign of the weight that depends on the configuration.
    This weight is equal to the `CGAL::Generalized_weights::Tangent_weight_2`.
    This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight_2`.

    \cgalFigureBegin{mean_value_weight, mean_value.svg}
      Notation used for the mean value weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Mean_value_weight_2 {

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
    Mean_value_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D mean value weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vm, vj, vp);
    }

    /*!
      \brief computes 2D mean value weight.
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

      const auto dot_product_2 =
        m_traits.compute_scalar_product_2_object();
      const FT Dm = dot_product_2(sm, sj);
      const FT Dj = dot_product_2(sj, sp);
      const FT Dq = dot_product_2(sm, sp);

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT Am = area_2(vm, vj, query);
      const FT Aj = area_2(vj, vp, query);
      const FT Bj = area_2(vm, vp, query);
      const FT sign = sign_of_weight(Am, Aj, Bj);

      return weight(
        rm, rj, rp, Dm, Dj, Dq, sign);
    }

    const FT weight_3(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      Point_2 pq, pm, pj, pp;
      internal::flatten(
        m_traits, query, vm, vj, vp,
        pq, pm, pj, pp);
      return weight_2(pq, pm, pj, pp);
    }

    const FT weight(
      const FT rm, const FT rj, const FT rp,
      const FT Dm, const FT Dj, const FT Dq,
      const FT sign) const {

      const FT Pm = rm * rj + Dm;
      const FT Pj = rj * rp + Dj;

      FT w = FT(0);
      CGAL_assertion(Pm != FT(0) && Pj != FT(0));
      const FT prod = Pm * Pj;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = FT(2) * (rm * rp - Dq) * inv;
        CGAL_assertion(w >= FT(0));
        if (w < FT(0)) w = CGAL::abs(w);
        w = static_cast<FT>(
          CGAL::sqrt(CGAL::to_double(w)));
      }
      w *= FT(2); w *= sign;
      return w;
    }

    FT sign_of_weight(
      const FT& Am,
      const FT& Aj,
      const FT& Bj) const {

      if (Am > FT(0) && Aj > FT(0) && Bj <= FT(0)) return  FT(1);
      if (Am < FT(0) && Aj < FT(0) && Bj >= FT(0)) return -FT(1);
      if (Bj > FT(0)) return  FT(1);
      if (Bj < FT(0)) return -FT(1);
      return FT(0);
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H
