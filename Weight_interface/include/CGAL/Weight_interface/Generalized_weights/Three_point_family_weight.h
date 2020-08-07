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

#ifndef CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHT_H
#define CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief Three point family weight.

    The full weight is computed as

    \f$w = \frac{r_p^a A_m - r^a B + r_m^a A}{A_m A}\f$

    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 0\f$ this weight is equal to the
    `CGAL::Generalized_weights::Wachspress_weight` and
    `CGAL::Generalized_weights::Authalic_weight`.

    For \f$a = 1\f$ this weight is equal to the
    `CGAL::Generalized_weights::Mean_value_weight` and
    `CGAL::Generalized_weights::Tangent_weight`.

    For \f$a = 2\f$ this weight is equal to the
    `CGAL::Generalized_weights::Discrete_harmonic_weight` and
    `CGAL::Generalized_weights::Cotangent_weight`.

    \cgalFigureBegin{three_point_family_weight, three_point_family.svg}
      Notation used for the three point family weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Three_point_family_weight {

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

      \param a
      the power parameter.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.
    */
    Three_point_family_weight(
      const FT a = FT(1), // default is for mean value coordinates
      const GeomTraits traits = GeomTraits()) :
    m_p(a), m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the three point family weight.
    */
    const FT operator()(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      return weight_2(q, t, r, p);
    }

    /*!
      \brief computes the three point family weight.
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
    const FT m_p;
    const GeomTraits m_traits;

    const FT weight_2(
      const Point_2& q,
      const Point_2& t,
      const Point_2& r,
      const Point_2& p) const {

      const FT d1 = internal::distance_2(m_traits, q, t);
      const FT d2 = internal::distance_2(m_traits, q, r);
      const FT d3 = internal::distance_2(m_traits, q, p);

      const FT A1 = internal::area_2(m_traits, r, q, t);
      const FT A2 = internal::area_2(m_traits, p, q, r);
      const FT B  = internal::area_2(m_traits, p, q, t);

      return weight(
        d1, d2, d3, A1, A2, B);
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
      const FT A1, const FT A2, const FT B) const {

      FT w = FT(0);
      CGAL_assertion(A1 != FT(0) && A2 != FT(0));
      const FT prod = A1 * A2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        FT a1 = r1;
        FT b  = r2;
        FT a2 = r3;
        if (m_p != FT(1)) {
          a1 = internal::power(m_traits, r1, m_p);
          b  = internal::power(m_traits, r2, m_p);
          a2 = internal::power(m_traits, r3, m_p);
        }
        w = (a1 * A1 - b * B + a2 * A2) * inv;
      }
      return w;
    }
  };

  template<typename Point_2>
  decltype(auto) three_point_family_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    Three_point_family_weight<Traits> family;
    return family(q, t, r, p);
  }

  template<typename Point_3>
  decltype(auto) three_point_family_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    Three_point_family_weight<Traits> family;
    return family(q, t, r, p);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_THREE_POINT_FAMILY_WEIGHT_H
