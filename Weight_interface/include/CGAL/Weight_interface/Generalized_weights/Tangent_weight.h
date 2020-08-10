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

#ifndef CGAL_GENERALIZED_TANGENT_WEIGHT_H
#define CGAL_GENERALIZED_TANGENT_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief Tangent weight.

    The full weight is computed as

    \f$w = 2 \frac{t_m + t}{r}\f$, where \f$t_m = \frac{A_m}{r r_m + D_m}\f$ and
    \f$t = \frac{A}{r r_p + D_p}\f$

    and the half weight as

    \f$h = 2 \frac{t}{r}\f$

    with notations shown in the figure below and dot products

    \f$D_m = (v_j - q) \cdot (v_m - q)\f$ and
    \f$D_p = (v_j - q) \cdot (v_p - q)\f$.

    This weight is equal to the `CGAL::Generalized_weights::Mean_value_weight`.
    This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight`.

    \cgalFigureBegin{tangent_weight, tangent.svg}
      Notation used for the tangent weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Tangent_weight {

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
    Tangent_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief to be added
    */
    const FT distance(
      const Point_2& p,
      const Point_2& q) const {

      return distance_2(p, q);
    }

    /*!
      \brief to be added
    */
    const FT area(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return area_2(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT scalar_product(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return scalar_product_2(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT distance(
      const Point_3& p,
      const Point_3& q) const {

      return distance_3(p, q);
    }

    /*!
      \brief to be added
    */
    const FT area(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return area_3(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT scalar_product(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return scalar_product_3(p, q, r);
    }

    /*!
      \brief to be added
    */
    const FT tangent(
      const FT d, const FT l, const FT A, const FT D) const {

      return tangent_half_angle(d, l, A, D);
    }

    /*!
      \brief computes the half of the tangent weight.
    */
    const FT operator()(
      const FT d, const FT tan) const {

      return half_weight(d, tan);
    }

    /*!
      \brief computes the half of the tangent weight.
    */
    const FT operator()(
      const FT d, const FT l, const FT A, const FT D) const {

      const FT tan = tangent_half_angle(d, l, A, D);
      return operator()(d, tan);
    }

    /*!
      \brief computes the tangent weight.
    */
    const FT operator()(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      return weight_2(q, t, r, p);
    }

    /*!
      \brief computes the tangent weight.
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

    const FT distance_2(
      const Point_2& p,
      const Point_2& q) const {

      return internal::distance_2(m_traits, p, q);
    }

    const FT area_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return internal::area_2(m_traits, p, q, r);
    }

    const FT scalar_product_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      const auto s1 = r - q;
      const auto s2 = p - q;

      const auto scalar_product_2 =
        m_traits.compute_scalar_product_2_object();
      return scalar_product_2(s1, s2);
    }

    const FT distance_3(
      const Point_3& p,
      const Point_3& q) const {

      return internal::distance_3(m_traits, p, q);
    }

    const FT area_3(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return internal::positive_area_3(m_traits, p, q, r);
    }

    const FT scalar_product_3(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      const auto s1 = r - q;
      const auto s2 = p - q;

      const auto scalar_product_3 =
        m_traits.compute_scalar_product_3_object();
      return scalar_product_3(s1, s2);
    }

    const FT tangent_half_angle(
      const FT r, const FT d,
      const FT A, const FT D) const {

      FT t = FT(0);
      const FT P = r * d + D;
      CGAL_assertion(P != FT(0));
      if (P != FT(0)) {
        const FT inv = FT(2) / P;
        t = A * inv;
      }
      return t;
    }

    const FT half_weight(
      const FT r, const FT t) const {

      FT w = FT(0);
      CGAL_assertion(r != FT(0));
      if (r != FT(0)) {
        const FT inv = FT(2) / r;
        w = t * inv;
      }
      return w;
    }

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

      const FT A1 = internal::area_2(m_traits, r, q, t);
      const FT A2 = internal::area_2(m_traits, p, q, r);

      const auto dot_product_2 =
        m_traits.compute_scalar_product_2_object();
      const FT D1 = dot_product_2(s1, s2);
      const FT D2 = dot_product_2(s2, s3);

      return weight(
        l1, l2, l3, A1, A2, D1, D2);
    }

    const FT weight_3(
      const Point_3& q,
      const Point_3& t,
      const Point_3& r,
      const Point_3& p) const {

      const auto s1 = t - q;
      const auto s2 = r - q;
      const auto s3 = p - q;

      const FT l1 = internal::length_3(m_traits, s1);
      const FT l2 = internal::length_3(m_traits, s2);
      const FT l3 = internal::length_3(m_traits, s3);

      const FT A1 = internal::positive_area_3(m_traits, r, q, t);
      const FT A2 = internal::positive_area_3(m_traits, p, q, r);

      const auto dot_product_3 =
        m_traits.compute_scalar_product_3_object();
      const FT D1 = dot_product_3(s1, s2);
      const FT D2 = dot_product_3(s2, s3);

      return weight(
        l1, l2, l3, A1, A2, D1, D2);
    }

    const FT weight(
      const FT d1, const FT r, const FT d2,
      const FT A1, const FT A2,
      const FT D1, const FT D2) const {

      const FT P1 = d1 * r + D1;
      const FT P2 = d2 * r + D2;

      FT w = FT(0);
      CGAL_assertion(P1 != FT(0) && P2 != FT(0));
      if (P1 != FT(0) && P2 != FT(0)) {
        const FT inv1 = FT(2) / P1;
        const FT inv2 = FT(2) / P2;
        const FT t1 = A1 * inv1;
        const FT t2 = A2 * inv2;
        CGAL_assertion(r != FT(0));
        if (r != FT(0)) {
          const FT inv = FT(2) / r;
          w = (t1 + t2) * inv;
        }
      }
      return w;
    }
  };

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the tangent weight for 2D points.

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
  decltype(auto) tangent_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    const Tangent_weight<Traits> tangent;
    return tangent(q, t, r, p);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the tangent weight for 3D points.

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
  decltype(auto) tangent_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    const Tangent_weight<Traits> tangent;
    return tangent(q, t, r, p);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_TANGENT_WEIGHT_H
