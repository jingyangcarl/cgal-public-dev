// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_2_H
#define CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos.
// A general construction of barycentric coordinates over convex polygons.
// Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  class Wachspress_weights_2 {

  public:
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    using Area_2 = typename GeomTraits::Compute_area_2;

    using FT = typename GeomTraits::FT;
    using Point_2 = typename GeomTraits::Point_2;

    Wachspress_weights_2(
      const Polygon& polygon,
      const GeomTraits traits,
      const VertexMap vertex_map) :
    m_polygon(polygon),
    m_traits(traits),
    m_vertex_map(vertex_map),
    m_area_2(m_traits.compute_area_2_object()) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, vertex_map));
      CGAL_precondition(
        internal::polygon_type_2(polygon, traits, vertex_map) ==
        internal::Polygon_type::STRICTLY_CONVEX);
      resize();
    }

    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      return optimal_weights(
        query, weights, normalize);
    }

  private:

    // Fields.
    const Polygon& m_polygon;
    const GeomTraits m_traits;
    const VertexMap m_vertex_map;
    const Area_2 m_area_2;

    std::vector<FT> A;
    std::vector<FT> C;
    std::vector<FT> w;

    // Functions.
    void resize() {
      A.resize(m_polygon.size());
      C.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator optimal_weights(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute areas A and C following the area notation from [1].
      // Split the loop to make this computation faster.
      const auto& p1 = get(m_vertex_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_vertex_map, *(m_polygon.begin() + 1));
      const auto& pn = get(m_vertex_map, *(m_polygon.begin() + (n - 1)));

      A[0] = m_area_2(p1, p2, query);
      C[0] = m_area_2(pn, p1, p2);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi0 = get(m_vertex_map, *(m_polygon.begin() + (i - 1)));
        const auto& pi1 = get(m_vertex_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_vertex_map, *(m_polygon.begin() + (i + 1)));

        A[i] = m_area_2(pi1, pi2, query);
        C[i] = m_area_2(pi0, pi1, pi2);
      }

      const auto& pm = get(m_vertex_map, *(m_polygon.begin() + (n - 2)));
      A[n - 1] = m_area_2(pn, p1, query);
      C[n - 1] = m_area_2(pm, pn, p1);

      // Compute unnormalized weights following the formula (28) from [1].
      CGAL_assertion(A[n - 1] != FT(0) && A[0] != FT(0));
      w[0] = C[0] / (A[n - 1] * A[0]);

      for (std::size_t i = 1; i < n - 1; ++i) {
        CGAL_assertion(A[i - 1] != FT(0) && A[i] != FT(0));
        w[i] = C[i] / (A[i - 1] * A[i]);
      }

      CGAL_assertion(A[n - 2] != FT(0) && A[n - 1] != FT(0));
      w[n - 1] = C[n - 1] / (A[n - 2] * A[n - 1]);

      // Normalize if necessary.
      if (normalize)
        internal::normalize(w);

      // Return weights.
      for (std::size_t i = 0; i < n; ++i)
        *(weights++) = w[i];
      return weights;
    }
  };

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_2_H
