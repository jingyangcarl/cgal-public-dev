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

#ifndef CGAL_GENERALIZED_WEIGHTS_INTERNAL_TOOLS_H
#define CGAL_GENERALIZED_WEIGHTS_INTERNAL_TOOLS_H

// #include <CGAL/license/Weight_interface.h>

// Boost includes.
#include <CGAL/boost/graph/helpers.h>

// Internal includes.
#include <CGAL/Weight_interface/Weighting_regions.h>
#include <CGAL/Weight_interface/Generalized_weights.h>

namespace CGAL {
namespace Generalized_weights {
namespace internal {

// Computes a secure cotanget between two 3D vectors.
template<typename Point_3>
decltype(auto) cotangent_3_secure(
  const Point_3& p,
  const Point_3& q,
  const Point_3& r) {

  using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
  using FT = typename GeomTraits::FT;
  using Get_sqrt = internal::Get_sqrt<GeomTraits>;
  const GeomTraits traits;
  const auto sqrt = Get_sqrt::sqrt_object(traits);

  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();
  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();
  const auto construct_vector_3 =
    traits.construct_vector_3_object();

  const auto v1 = construct_vector_3(q, r);
  const auto v2 = construct_vector_3(q, p);

  const FT dot = dot_product_3(v1, v2);
  const FT length_v1 = length_3(traits, v1);
  const FT length_v2 = length_3(traits, v2);

  const FT lb = -FT(999) / FT(1000), ub = FT(999) / FT(1000);
  FT cosine = dot / length_v1 / length_v2;
  cosine = (cosine < lb) ? lb : cosine;
  cosine = (cosine > ub) ? ub : cosine;
  const FT sine = sqrt(FT(1) - cosine * cosine);

  CGAL_assertion(sine != FT(0));
  if (sine != FT(0)) return cosine / sine;
  return FT(0); // undefined
}

template<typename GeomTraits>
class Tangent_weight_wrapper {

  using FT = typename GeomTraits::FT;
  FT m_d_r, m_d_p, m_w_base;

public:
  template<typename Point>
  Tangent_weight_wrapper(
    const Point& p,
    const Point& q,
    const Point& r) {

    m_d_r = CGAL::Generalized_weights::utils::distance_3(q, r);
    CGAL_assertion(m_d_r != FT(0)); // two points are identical!
    m_d_p = CGAL::Generalized_weights::utils::distance_3(q, p);
    CGAL_assertion(m_d_p != FT(0)); // two points are identical!
    const FT area = CGAL::Generalized_weights::utils::area_3(p, q, r);
    CGAL_assertion(area != FT(0));  // three points are identical!
    const FT scalar = CGAL::Generalized_weights::utils::scalar_product_3(p, q, r);

    m_w_base = -CGAL::Generalized_weights::utils::
      tangent_half_angle(m_d_r, m_d_p, area, scalar);
  }

  const FT get_w_r() const {
    return CGAL::Generalized_weights::utils::
      half_tangent_weight(m_w_base, m_d_r) / FT(2);
  }

  const FT get_w_p() const {
    return CGAL::Generalized_weights::utils::
      half_tangent_weight(m_w_base, m_d_p) / FT(2);
  }
};

template<
typename GeomTraits,
typename PolygonMesh>
class Cotangent_weight_wrapper {

  using FT = typename GeomTraits::FT;
  bool m_use_secure_version;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Cotangent_weight_wrapper(
    const bool use_secure_version = false) :
  m_use_secure_version(use_secure_version)
  { }

  template<class VertexPointMap>
  FT operator()(
    const PolygonMesh& pmesh,
    const halfedge_descriptor he,
    const VertexPointMap& pmap) const {

    const auto v0 = target(he, pmesh);
    const auto v1 = source(he, pmesh);

    const auto& p0 = get(pmap, v0);
    const auto& p1 = get(pmap, v1);

    FT weight = FT(0);
    if (is_border_edge(he, pmesh)) {
      const auto he_cw = opposite(next(he, pmesh), pmesh);
      auto v2 = source(he_cw, pmesh);

      if (is_border_edge(he_cw, pmesh)) {
        const auto he_ccw = prev(opposite(he, pmesh), pmesh);
        v2 = source(he_ccw, pmesh);

        const auto& p2 = get(pmap, v2);
        if (m_use_secure_version)
          weight = internal::cotangent_3_secure(p1, p2, p0);
        else
          weight = CGAL::Generalized_weights::utils::cotangent_3(p1, p2, p0);
        weight = (CGAL::max)(FT(0), weight);
        weight /= FT(2);
      } else {
        const auto& p2 = get(pmap, v2);
        if (m_use_secure_version)
          weight = internal::cotangent_3_secure(p0, p2, p1);
        else
          weight = CGAL::Generalized_weights::utils::cotangent_3(p0, p2, p1);
        weight = (CGAL::max)(FT(0), weight);
        weight /= FT(2);
      }

    } else {
      const auto he_cw = opposite(next(he, pmesh), pmesh);
      const auto v2 = source(he_cw, pmesh);
      const auto he_ccw = prev(opposite(he, pmesh), pmesh);
      const auto v3 = source(he_ccw, pmesh);

      const auto& p2 = get(pmap, v2);
      const auto& p3 = get(pmap, v3);
      FT cot_beta = FT(0), cot_gamma = FT(0);

      if (m_use_secure_version)
        cot_beta = internal::cotangent_3_secure(p0, p2, p1);
      else
        cot_beta = CGAL::Generalized_weights::utils::cotangent_3(p0, p2, p1);

      if (m_use_secure_version)
        cot_gamma = internal::cotangent_3_secure(p1, p3, p0);
      else
        cot_gamma = CGAL::Generalized_weights::utils::cotangent_3(p1, p3, p0);

      cot_beta  = (CGAL::max)(FT(0), cot_beta);  cot_beta  /= FT(2);
      cot_gamma = (CGAL::max)(FT(0), cot_gamma); cot_gamma /= FT(2);
      weight = cot_beta + cot_gamma;
    }
    return weight;
  }
};

template<
typename PolygonMesh,
typename VertexPointMap = typename boost::property_map<PolygonMesh, vertex_point_t>::type>
class Cotangent_weight_wrapper_with_voronoi_secure {

  using GeomTraits = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
  using FT = typename GeomTraits::FT;

  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  Cotangent_weight_wrapper_with_voronoi_secure(
    const PolygonMesh& pmesh,
    const VertexPointMap pmap) :
  m_pmesh(pmesh),
  m_pmap(pmap)
  { }

  PolygonMesh& pmesh() {
    return m_pmesh;
  }

  const FT w_i(const vertex_descriptor v_i) const {
    return FT(1) / (FT(2) * voronoi(v_i));
  }

  const FT w_ij(const halfedge_descriptor he) const {
    return cotangent_secure(he);
  }

private:
  const FT cotangent_secure(halfedge_descriptor he) const {

    const auto v0 = target(he, m_pmesh);
    const auto v1 = source(he, m_pmesh);

    const auto& p0 = get(m_pmap, v0);
    const auto& p1 = get(m_pmap, v1);

    FT weight = FT(0);
    if (is_border_edge(he, m_pmesh)) {
      const auto he_cw = opposite(next(he, m_pmesh), m_pmesh);
      auto v2 = source(he_cw, m_pmesh);

      if (is_border_edge(he_cw, m_pmesh)) {
        const auto he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
        v2 = source(he_ccw, m_pmesh);

        const auto& p2 = get(m_pmap, v2);
        weight = internal::cotangent_3_secure(p1, p2, p0);
      } else {
        const auto& p2 = get(m_pmap, v2);
        weight = internal::cotangent_3_secure(p0, p2, p1);
      }

    } else {
      const auto he_cw = opposite(next(he, m_pmesh), m_pmesh);
      const auto v2 = source(he_cw, m_pmesh);
      const auto he_ccw = prev(opposite(he, m_pmesh), m_pmesh);
      const auto v3 = source(he_ccw, m_pmesh);

      const auto& p2 = get(m_pmap, v2);
      const auto& p3 = get(m_pmap, v3);

      const FT cot_beta  = internal::cotangent_3_secure(p0, p2, p1);
      const FT cot_gamma = internal::cotangent_3_secure(p1, p3, p0);
      weight = cot_beta + cot_gamma;
    }
    return weight;
  }

  const FT voronoi(const vertex_descriptor v0) const {

    FT voronoi_area = FT(0);
    CGAL_assertion(CGAL::is_triangle_mesh(m_pmesh));
    for (const auto he :
      halfedges_around_target(halfedge(v0, m_pmesh), m_pmesh)) {
      if (is_border(he, m_pmesh)) continue;
      CGAL_assertion(v0 == target(he, m_pmesh));

      const auto v1 = source(he, m_pmesh);
      const auto v2 = target(next(he, m_pmesh), m_pmesh);

      const auto& p0 = get(m_pmap, v0);
      const auto& p1 = get(m_pmap, v1);
      const auto& p2 = get(m_pmap, v2);

      voronoi_area += CGAL::Generalized_weights::
        mixed_voronoi_area_3(p1, p0, p2);
    }
    CGAL_assertion(voronoi_area != FT(0));
    return voronoi_area;
  }
};

template<
typename GeomTraits,
typename PolygonMesh,
typename VertexPointMap>
class Edge_cotangent_weight_wrapper {

  using FT = typename GeomTraits::FT;
  const PolygonMesh& m_pmesh;
  const VertexPointMap m_pmap;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  PM_edge_cotangent_weight(
    PolygonMesh& pmesh, VertexPointMap pmap) :
  m_pmesh(pmesh),
  m_pmap(pmap)
  { }

  PolygonMesh& pmesh() {
    return m_pmesh;
  }

  const FT operator()(const halfedge_descriptor he) const {

    FT weight = FT(0);
    if (is_border_edge(he, m_pmesh)) {
      const auto h1 = next(he, m_pmesh);

      const auto v0 = target(he, m_pmesh);
      const auto v1 = source(he, m_pmesh);
      const auto v2 = target(h1, m_pmesh);

      const auto& p0 = get(m_pmap, v0);
      const auto& p1 = get(m_pmap, v1);
      const auto& p2 = get(m_pmap, v2);

      weight = CGAL::Generalized_weights::utils::
        cotangent_3(p0, p2, p1);

    } else {
      const auto h1 = next(he, m_pmesh);
      const auto h2 = prev(opposite(he, m_pmesh), m_pmesh);

      const auto v0 = target(he, m_pmesh);
      const auto v1 = source(he, m_pmesh);
      const auto v2 = target(h1, m_pmesh);
      const auto v3 = source(h2, m_pmesh);

      const auto& p0 = get(m_pmap, v0);
      const auto& p1 = get(m_pmap, v1);
      const auto& p2 = get(m_pmap, v2);
      const auto& p3 = get(m_pmap, v3);

      weight = CGAL::Generalized_weights::
        cotangent_weight_3(p0, p2, p1, p3) / FT(2);
    }
    return weight;
  }
};

template<
typename GeomTraits,
typename PolygonMesh>
class Single_cotangent_weight_wrapper {

  using FT = typename GeomTraits::FT;

public:
  using halfedge_descriptor = typename boost::graph_traits<PolygonMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;

  template<class VertexPointMap>
  FT operator()(
    const PolygonMesh& pmesh,
    const halfedge_descriptor he,
    const VertexPointMap& pmap) const {

    if (is_border(he, pmesh)) return FT(0);

    const vertex_descriptor v0 = target(he, pmesh);
    const vertex_descriptor v1 = source(he, pmesh);
    const vertex_descriptor v2 = target(next(he, pmesh), pmesh);

    const auto& p0 = get(pmap, v0);
    const auto& p1 = get(pmap, v1);
    const auto& p2 = get(pmap, v2);

    return CGAL::Generalized_weights::utils::cotangent_3(p0, p2, p1);
  }
};

} // namespace internal
} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_INTERNAL_TOOLS_H
