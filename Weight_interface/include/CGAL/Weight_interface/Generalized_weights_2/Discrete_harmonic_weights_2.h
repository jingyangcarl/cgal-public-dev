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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHTS_2_H
#define CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHTS_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/Generalized_weights_2/enum_2.h>
#include <CGAL/Weight_interface/Generalized_weights_2/internal/utils_2.h>

// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos.
// A general construction of barycentric coordinates over convex polygons.
// Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2D

    \brief 2D discrete harmonic weights.

    This class implements 2D discrete harmonic weights ( \cite cgal:bc:fhk-gcbcocp-06,
    \cite cgal:pp-cdmsc-93, \cite cgal:bc:eddhls-maam-95 ).

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits
    is a model of `AnalyticTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \cgalModels `AnalyticWeights_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Discrete_harmonic_weights_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    using Area_2 = typename GeomTraits::Compute_area_2;
    using Squared_distance_2 = typename GeomTraits::Compute_squared_distance_2;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      This class implements the behavior of discrete harmonic weights
      for 2D query points.

      \param polygon
      An instance of `Polygon` with the vertices of a strictly convex polygon.

      \param policy
      One of the `CGAL::Generalized_weights::Computation_policy_2`.
      The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon`
      to `Point_2`. The default is the identity property map.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
      \pre `polygon is strictly convex`
    */
    Discrete_harmonic_weights_2(
      const Polygon& polygon,
      const Computation_policy_2 policy
      = Computation_policy_2::DEFAULT,
      const GeomTraits traits = GeomTraits(),
      const VertexMap vertex_map = VertexMap()) :
    m_polygon(polygon),
    m_computation_policy(policy),
    m_traits(traits),
    m_vertex_map(vertex_map),
    m_area_2(m_traits.compute_area_2_object()),
    m_squared_distance_2(m_traits.compute_squared_distance_2_object()) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, vertex_map));
      CGAL_precondition(
        internal::polygon_type_2(polygon, traits, vertex_map) ==
        internal::Polygon_type::STRICTLY_CONVEX);
      resize();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D discrete harmonic weights.

      This function fills `weights` with 2D discrete harmonic weights computed at the `query`
      point with respect to the vertices of the input polygon. If `query` belongs to
      the polygon boundary, the returned weights are normalized.

      The number of returned weights equals to the number of polygon vertices.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param weights
      An output iterator that stores the computed weights.

      \return an output iterator.
    */
    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator weights) {

      const bool normalize = false;
      return compute(query, weights, normalize);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      return compute(query, weights, normalize);
    }
    /// \endcond

  private:

    // Fields.
    const Polygon& m_polygon;
    const Computation_policy_2 m_computation_policy;
    const GeomTraits m_traits;
    const VertexMap m_vertex_map;
    const Area_2 m_area_2;
    const Squared_distance_2 m_squared_distance_2;

    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> B;
    std::vector<FT> w;

    // Functions.
    void resize() {
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      B.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator compute(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      switch (m_computation_policy) {

        case Computation_policy_2::OPTIMAL: {
          return optimal_weights(query, weights, normalize);
        }

        case Computation_policy_2::OPTIMAL_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          if (edge_case == internal::Edge_case::EXTERIOR)
            std::cerr << std::endl <<
            "WARNING: query does not belong to the polygon!" << std::endl;
          return optimal_weights(query, weights, normalize);
        }

        default: {
          internal::get_default(m_polygon.size(), weights);
          return weights;
        }
      }
      return weights;
    }

    template<typename OutputIterator>
    internal::Edge_case verify(
      const Point_2& query,
      OutputIterator weights) const {

      const auto result = internal::locate_wrt_polygon_2(
        m_polygon, query, m_traits, m_vertex_map);
      if (!result)
        return internal::Edge_case::EXTERIOR;

      const auto location = (*result).first;
      const std::size_t index = (*result).second;
      if (location == internal::Query_point_location::ON_UNBOUNDED_SIDE)
        return internal::Edge_case::EXTERIOR;

      if (
        location == internal::Query_point_location::ON_VERTEX ||
        location == internal::Query_point_location::ON_EDGE ) {
        internal::boundary_coordinates_2(
          m_polygon, query, location, index, weights, m_traits, m_vertex_map);
        return internal::Edge_case::BOUNDARY;
      }
      return internal::Edge_case::INTERIOR;
    }

    template<typename OutputIterator>
    OutputIterator optimal_weights(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute areas A, B, and distances r following the notation from [1].
      // Split the loop to make this computation faster.
      const auto& p1 = get(m_vertex_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_vertex_map, *(m_polygon.begin() + 1));
      const auto& pn = get(m_vertex_map, *(m_polygon.begin() + (n - 1)));

      r[0] = m_squared_distance_2(p1, query);
      A[0] = m_area_2(p1, p2, query);
      B[0] = m_area_2(pn, p2, query);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi0 = get(m_vertex_map, *(m_polygon.begin() + (i - 1)));
        const auto& pi1 = get(m_vertex_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_vertex_map, *(m_polygon.begin() + (i + 1)));

        r[i] = m_squared_distance_2(pi1, query);
        A[i] = m_area_2(pi1, pi2, query);
        B[i] = m_area_2(pi0, pi2, query);
      }

      const auto& pm = get(m_vertex_map, *(m_polygon.begin() + (n - 2)));
      r[n - 1] = m_squared_distance_2(pn, query);
      A[n - 1] = m_area_2(pn, p1, query);
      B[n - 1] = m_area_2(pm, p1, query);

      // Compute unnormalized weights following the formula (25) with p = 2 from [1].
      CGAL_precondition(A[n - 1] != FT(0) && A[0] != FT(0));
      w[0] = (r[1] * A[n - 1] - r[0] * B[0] + r[n - 1] * A[0]) / (A[n - 1] * A[0]);

      for (std::size_t i = 1; i < n - 1; ++i) {
        CGAL_precondition(A[i - 1] != FT(0) && A[i] != FT(0));
        w[i] = (r[i + 1] * A[i - 1] - r[i] * B[i] + r[i - 1] * A[i]) / (A[i - 1] * A[i]);
      }

      CGAL_precondition(A[n - 2] != FT(0) && A[n - 1] != FT(0));
      w[n - 1] = (r[0] * A[n - 2] - r[n - 1] * B[n - 1] + r[n - 2] * A[n - 1]) / (A[n - 2] * A[n - 1]);

      // Normalize if necessary.
      if (normalize)
        internal::normalize(w);

      // Return weights.
      for (std::size_t i = 0; i < n; ++i)
        *(weights++) = w[i];

      // Return weights.
      return weights;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_DISCRETE_HARMONIC_WEIGHTS_2_H
