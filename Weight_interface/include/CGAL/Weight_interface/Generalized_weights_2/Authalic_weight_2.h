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

#ifndef CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H
#define CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief 2D authalic weight.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`, `HalfWeight_2`
  */
  template<typename GeomTraits>
  class Authalic_weight_2 {

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
    Authalic_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the half value of a 2D authalic weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      return half_weight_2(query, vj, vp);
    }

    /*!
      \brief computes the half value of a 2D authalic weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      return half_weight_3(query, vj, vp);
    }

    /*!
      \brief computes 2D authalic weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vm, vj, vp);
    }

    /*!
      \brief computes 2D authalic weight.
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

      const FT cot_angle =
        internal::cotangent_2(m_traits, query, vj, vp);
      const auto squared_distance_2 =
        m_traits.compute_squared_distance_2_object();
      const FT rj2 = squared_distance_2(query, vj);
      return half_weight(cot_angle, rj2);
    }

    const FT half_weight_3(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      const FT cot_angle =
        internal::cotangent_3(m_traits, query, vj, vp);
      const auto squared_distance_3 =
        m_traits.compute_squared_distance_3_object();
      const FT rj2 = squared_distance_3(query, vj);
      return half_weight(cot_angle, rj2);
    }

    const FT half_weight(
      const FT cot_angle, const FT rj2) const {

      FT w = FT(0);
      CGAL_assertion(rj2 != FT(0));
      if (rj2 != FT(0)) {
        const FT inv = FT(2) / rj2;
        w = cot_angle * inv;
      }
      return w;
    }

    const FT weight_2(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      const FT cot_gamma = internal::cotangent_2(m_traits, vm, vj, query);
      const FT cot_beta  = internal::cotangent_2(m_traits, query, vj, vp);

      const auto squared_distance_2 =
        m_traits.compute_squared_distance_2_object();
      const FT rj2 = squared_distance_2(query, vj);

      return weight(
        cot_gamma, cot_beta, rj2);
    }

    const FT weight_3(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      const FT cot_gamma = internal::cotangent_3(m_traits, vm, vj, query);
      const FT cot_beta  = internal::cotangent_3(m_traits, query, vj, vp);

      const auto squared_distance_3 =
        m_traits.compute_squared_distance_3_object();
      const FT rj2 = squared_distance_3(query, vj);

      return weight(
        cot_gamma, cot_beta, rj2);
    }

    const FT weight(
      const FT cot_gamma,
      const FT cot_beta,
      const FT rj2) const {

      FT w = FT(0);
      CGAL_assertion(rj2 != FT(0));
      if (rj2 != FT(0)) {
        const FT inv = FT(2) / rj2;
        w = (cot_gamma + cot_beta) * inv;
      }
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H
