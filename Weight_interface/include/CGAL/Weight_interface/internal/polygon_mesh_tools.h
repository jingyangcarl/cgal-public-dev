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

#ifndef CGAL_GENERALIZED_WEIGHTS_POLYGON_MESH_TOOLS_H
#define CGAL_GENERALIZED_WEIGHTS_POLYGON_MESH_TOOLS_H

// #include <CGAL/license/Weight_interface.h>

// Boost includes.
#include <CGAL/boost/graph/helpers.h>

// Internal includes.
#include <CGAL/Weight_interface/Generalized_weights/Tangent_weight.h>
#include <CGAL/Weight_interface/Generalized_weights/Cotangent_weight.h>

namespace CGAL {
namespace Generalized_weights {
namespace internal {

template<
typename GeomTraits,
typename PolygonMesh>
class PM_tangent_weight {

  using FT = typename GeomTraits::FT;
  using Tangent_weight = CGAL::Generalized_weights::Tangent_weight<GeomTraits>;
  const Tangent_weight m_tangent_weight;
  FT m_d_r, m_d_p, m_w_base;

public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;

  PM_tangent_weight() { }

  template<typename Point>
  PM_tangent_weight(
    const Point& p,
    const Point& q,
    const Point& r) {

    m_d_r = m_tangent_weight.distance(q, r);
    CGAL_assertion(m_d_r != FT(0)); // two points are identical!
    m_d_p = m_tangent_weight.distance(q, p);
    CGAL_assertion(m_d_p != FT(0)); // two points are identical!
    const FT area = m_tangent_weight.area(p, q, r);
    CGAL_assertion(area != FT(0));  // three points are identical!
    const FT scalar = m_tangent_weight.scalar_product(p, q, r);

    m_w_base = -m_tangent_weight.tangent(
      m_d_r, m_d_p, area, scalar);
  }

  const FT get_w_r() const {
    return m_tangent_weight(m_d_r, m_w_base) / FT(2);
  }

  const FT get_w_p() const {
    return m_tangent_weight(m_d_p, m_w_base) / FT(2);
  }
};

template<
typename GeomTraits,
typename PolygonMesh>
class PM_cotangent_weight {

  using FT = typename GeomTraits::FT;
  using Cotangent_weight = CGAL::Generalized_weights::Cotangent_weight<GeomTraits>;
  const Cotangent_weight m_cotangent_weight;
  bool m_use_secure_version;
  GeomTraits m_traits;

public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;

  PM_cotangent_weight(
    const bool use_secure_version = false,
    GeomTraits traits = GeomTraits()) :
  m_use_secure_version(use_secure_version),
  m_traits(traits)
  { }

  template<class VertexPointMap>
  FT operator()(
    halfedge_descriptor he,
    PolygonMesh& pmesh,
    const VertexPointMap& ppmap) const {

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const auto& p0 = get(ppmap, v0);
    const auto& p1 = get(ppmap, v1);

    FT weight = FT(0);
    if (is_border_edge(he, pmesh)) {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh), pmesh);
      vertex_descriptor v2 = source(he_cw, pmesh);

      if (is_border_edge(he_cw, pmesh)) {
        const halfedge_descriptor he_ccw = prev(opposite(he, pmesh), pmesh);
        v2 = source(he_ccw, pmesh);

        const auto& p2 = get(ppmap, v2);
        if (m_use_secure_version)
          weight = internal::cotangent_3_secure(m_traits, p1, p2, p0);
        else
          weight = m_cotangent_weight.cotangent(p1, p2, p0);
        weight = (CGAL::max)(FT(0), weight);
        weight /= 2.0;
      } else {
        const auto& p2 = get(ppmap, v2);
        if (m_use_secure_version)
          weight = internal::cotangent_3_secure(m_traits, p0, p2, p1);
        else
          weight = m_cotangent_weight.cotangent(p0, p2, p1);
        weight = (CGAL::max)(FT(0), weight);
        weight /= 2.0;
      }

    } else {
      const halfedge_descriptor he_cw = opposite(next(he, pmesh), pmesh);
      const vertex_descriptor v2 = source(he_cw, pmesh);
      const halfedge_descriptor he_ccw = prev(opposite(he, pmesh), pmesh);
      const vertex_descriptor v3 = source(he_ccw, pmesh);

      const auto& p0 = get(ppmap, v0);
      const auto& p1 = get(ppmap, v1);
      const auto& p2 = get(ppmap, v2);
      const auto& p3 = get(ppmap, v3);

      FT cot_beta = FT(0), cot_gamma = FT(0);

      if (m_use_secure_version)
        cot_beta = internal::cotangent_3_secure(m_traits, p1, p3, p0);
      else
        cot_beta = m_cotangent_weight.cotangent(p1, p3, p0);

      if (m_use_secure_version)
        cot_gamma = internal::cotangent_3_secure(m_traits, p0, p2, p1);
      else
        cot_gamma = m_cotangent_weight.cotangent(p0, p2, p1);

      cot_beta  = (CGAL::max)(FT(0), cot_beta);  cot_beta  /= 2.0;
      cot_gamma = (CGAL::max)(FT(0), cot_gamma); cot_gamma /= 2.0;
      weight = cot_beta + cot_gamma;
    }
    return weight;
  }
};

template<
typename GeomTraits,
typename PolygonMesh>
class PM_single_cotangent_weight {

  using FT = typename GeomTraits::FT;
  using Cotangent_weight = CGAL::Generalized_weights::Cotangent_weight<GeomTraits>;
  const Cotangent_weight m_cotangent_weight;

public:
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;

  template<class VertexPointMap>
  FT operator()(
    halfedge_descriptor he,
    PolygonMesh& pmesh,
    const VertexPointMap& ppmap) const {

    if (is_border(he, pmesh)) { return FT(0); }

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const vertex_descriptor v2 = target(next(he, pmesh), pmesh);
    const auto& p0 = get(ppmap, v0);
    const auto& p1 = get(ppmap, v1);
    const auto& p2 = get(ppmap, v2);

    return m_cotangent_weight.cotangent(p0, p2, p1);
  }
};

} // namespace internal
} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_POLYGON_MESH_TOOLS_H
