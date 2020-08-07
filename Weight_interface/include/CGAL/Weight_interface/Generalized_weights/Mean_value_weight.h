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

#ifndef CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H
#define CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief Mean value weight.

    The full weight is computed as

    \f$w = \pm 2 \sqrt{\frac{2 (r_m r_p - D)}{(r r_m + D_m)(r r_p + D_p)}}\f$,

    with notations shown in the figure below and dot products

    \f$D_m = (v_j - q) \cdot (v_m - q)\f$,
    \f$D_p = (v_j - q) \cdot (v_p - q)\f$, and
    \f$D   = (v_m - q) \cdot (v_p - q)\f$.

    The \f$\pm\f$ sign is a sign of the weight that depends on the configuration.
    This weight is equal to the `CGAL::Generalized_weights::Tangent_weight`.
    This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight`.

    \cgalFigureBegin{mean_value_weight, mean_value.svg}
      Notation used for the mean value weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Mean_value_weight {

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
    Mean_value_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the mean value weight.
    */
    const FT operator()(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      return weight_2(q, t, r, p);
    }

    /*!
      \brief computes the mean value weight.
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

    const FT weight_2(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      const auto s1 = t - q;
      const auto s2 = r - q;
      const auto s3 = p - q;

      const FT l1 = internal::length_2(m_traits, s1);
      const FT l2 = internal::length_2(m_traits, s2);
      const FT l3 = internal::length_2(m_traits, s3);

      const auto dot_product_2 =
        m_traits.compute_scalar_product_2_object();
      const FT D1 = dot_product_2(s1, s2);
      const FT D2 = dot_product_2(s2, s3);
      const FT D  = dot_product_2(s1, s3);

      const FT A1 = internal::area_2(m_traits, r, q, t);
      const FT A2 = internal::area_2(m_traits, p, q, r);
      const FT B  = internal::area_2(m_traits, p, q, t);
      const FT sign = sign_of_weight(A1, A2, B);

      return weight(
        l1, l2, l3, D1, D2, D, sign);
    }

    const FT weight_3(
      const Point_3& q,
      const Point_3& t,
      const Point_3& r,
      const Point_3& p) const {

      Point_2 qf, tf, rf, pf;
      internal::flatten(
        m_traits,
        q,  t,  r,  p,
        qf, tf, rf, pf);
      return weight_2(qf, tf, rf, pf);
    }

    const FT weight(
      const FT r1, const FT r2, const FT r3,
      const FT D1, const FT D2, const FT D,
      const FT sign) const {

      const FT P1 = r1 * r2 + D1;
      const FT P2 = r2 * r3 + D2;

      FT w = FT(0);
      CGAL_assertion(P1 != FT(0) && P2 != FT(0));
      const FT prod = P1 * P2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = FT(2) * (r1 * r3 - D) * inv;
        CGAL_assertion(w >= FT(0));
        if (w < FT(0)) w = CGAL::abs(w);
        w = static_cast<FT>(
          CGAL::sqrt(CGAL::to_double(w)));
      }
      w *= FT(2); w *= sign;
      return w;
    }

    FT sign_of_weight(
      const FT& A1, const FT& A2, const FT& B) const {

      if (A1 > FT(0) && A2 > FT(0) && B <= FT(0)) return  FT(1);
      if (A1 < FT(0) && A2 < FT(0) && B >= FT(0)) return -FT(1);
      if (B  > FT(0)) return  FT(1);
      if (B  < FT(0)) return -FT(1);
      return FT(0);
    }
  };

  template<typename Point_2>
  decltype(auto) mean_value_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    Mean_value_weight<Traits> mean_value;
    return mean_value(q, t, r, p);
  }

  template<typename Point_3>
  decltype(auto) mean_value_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    Mean_value_weight<Traits> mean_value;
    return mean_value(q, t, r, p);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_H
