// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHT_2_H
#define CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Boost includes.
#include <CGAL/boost/graph/helpers.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2D

    \brief 2D discrete harmonic weight.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Discrete_harmonic_weight_2 {

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
    Discrete_harmonic_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D discrete harmonic weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vm, vj, vp);
    }

    /*!
      \brief computes 2D discrete harmonic weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      return weight_3(query, vm, vj, vp);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<
    typename PolygonMesh,
    typename VertexDescriptor,
    typename VertextAroundTargetCirculator>
    const FT operator()(
      const PolygonMesh& polygon_mesh,
      const VertexDescriptor vdi,
      const VertextAroundTargetCirculator vcj) const {

      const auto point_map = get(vertex_point, polygon_mesh);
      const Point_3& query = get(point_map, vdi);

      auto vcm = vcj; vcm--;
      auto vcp = vcj; vcp++;

      const Point_3& vm = get(point_map, vcm);
      const Point_3& vj = get(point_map, vcj);
      const Point_3& vp = get(point_map, vcp);

      return weight_3(query, vm, vj, vp);
    }
    /// \endcond

  private:
    const GeomTraits m_traits;

    const FT weight_2(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      const auto squared_distance_2 =
        m_traits.compute_squared_distance_2_object();
      const FT rm2 = squared_distance_2(query, vm);
      const FT rj2 = squared_distance_2(query, vj);
      const FT rp2 = squared_distance_2(query, vp);

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT Am = area_2(vm, vj, query);
      const FT Aj = area_2(vj, vp, query);
      const FT Bj = area_2(vm, vp, query);

      return weight(
        rm2, rj2, rp2, Am, Aj, Bj);
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
      const FT rm2, const FT rj2, const FT rp2,
      const FT Am, const FT Aj, const FT Bj) const {

      FT w = FT(0);
      CGAL_assertion(Am != FT(0) && Aj != FT(0));
      const FT prod = Am * Aj;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = (rp2 * Am - rj2 * Bj + rm2 * Aj) * inv;
      }
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHT_2_H
