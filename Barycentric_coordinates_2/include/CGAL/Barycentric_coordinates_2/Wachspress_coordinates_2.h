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

#ifndef CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_2_H
#define CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>
#include <CGAL/Barycentric_coordinates_2/internal/Generalized_weights_2/Wachspress_weights_2.h>

// [1] Reference: "M. S. Floater, K. Hormann, and G. Kos.
// A general construction of barycentric coordinates over convex polygons.
// Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefAnalytic

    \brief 2D Wachspress coordinates.

    This class implements 2D Wachspress coordinates ( \cite cgal:bc:fhk-gcbcocp-06,
    \cite cgal:bc:mlbd-gbcip-02, \cite cgal:bc:w-rfeb-75 ), which can be computed
    at any point inside a strictly convex polygon.

    Wachspress coordinates are well-defined and non-negative in the closure
    of a strictly convex polygon. The coordinates are computed analytically.

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Wachspress_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    using Area_2 = typename GeomTraits::Compute_area_2;

    using Wachspress_weights_2 =
      CGAL::Generalized_weights::Wachspress_weights_2<Polygon, GeomTraits, Vertex_map>;
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

      This class implements the behavior of Wachspress coordinates
      for 2D query points.

      \param polygon
      An instance of `Polygon` with the vertices of a strictly convex polygon.

      \param policy
      One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
      The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon`
      to `Point_2`. The default is the identity property map.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
      \pre `polygon is strictly convex`
    */
    Wachspress_coordinates_2(
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
    m_wachspress_weights_2(
      polygon,
      CGAL::Generalized_weights::Computation_policy_2::OPTIMAL,
      traits,
      vertex_map) {

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
      \brief computes 2D Wachspress weights.

      This function fills `weights` with 2D Wachspress weights computed at the `query`
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
    OutputIterator weights(
      const Point_2& query,
      OutputIterator weights) {

      const bool normalize = false;
      return compute(query, weights, normalize);
    }

    /*!
      \brief computes 2D Wachspress coordinates.

      This function fills `coordinates` with 2D Wachspress coordinates computed at the `query`
      point with respect to the vertices of the input polygon.

      The number of returned coordinates equals to the number of polygon vertices.

      \tparam OutputIterator
      is an output iterator whose value type is `FT`.

      \param query
      A query point.

      \param coordinates
      An output iterator that stores the computed coordinates.

      \return an output iterator.
    */
    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator coordinates) {

      const bool normalize = true;
      return compute(query, coordinates, normalize);
    }

    /// @}

  private:

    // Fields.
    const Polygon& m_polygon;
    const Computation_policy_2 m_computation_policy;
    const GeomTraits m_traits;
    const VertexMap m_vertex_map;
    const Area_2 m_area_2;

    Wachspress_weights_2 m_wachspress_weights_2;

    std::vector<FT> A;
    std::vector<FT> w;

    // Functions.
    void resize() {
      A.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator compute(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      switch (m_computation_policy) {

        case Computation_policy_2::PRECISE_COMPUTATION: {
          if (normalize) {
            return max_precision_weights(query, weights, normalize);
          } else {
            std::cerr << "WARNING: you can't use the precise version of unnormalized weights! ";
            std::cerr << "They are not valid weights!" << std::endl;
            internal::get_default(m_polygon.size(), weights);
            return weights;
          }
        }

        case Computation_policy_2::PRECISE_COMPUTATION_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          if (edge_case == internal::Edge_case::EXTERIOR)
            std::cerr << std::endl <<
            "WARNING: query does not belong to the polygon!" << std::endl;
          if (normalize) {
            return max_precision_weights(query, weights, normalize);
          } else {
            std::cerr << "WARNING: you can't use the precise version of unnormalized weights! ";
            std::cerr << "They are not valid weights!" << std::endl;
            internal::get_default(m_polygon.size(), weights);
            return weights;
          }
        }

        case Computation_policy_2::FAST_COMPUTATION: {
          return m_wachspress_weights_2(query, weights, normalize);
        }

        case Computation_policy_2::FAST_COMPUTATION_WITH_EDGE_CASES: {
          const auto edge_case = verify(query, weights);
          if (edge_case == internal::Edge_case::BOUNDARY)
            return weights;
          if (edge_case == internal::Edge_case::EXTERIOR)
            std::cerr << std::endl <<
            "WARNING: query does not belong to the polygon!" << std::endl;
          return m_wachspress_weights_2(query, weights, normalize);
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
    OutputIterator max_precision_weights(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute areas A following the area notation from [1].
      // Split the loop to make this computation faster.
      const auto& p1 = get(m_vertex_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_vertex_map, *(m_polygon.begin() + 1));
      A[0] = m_area_2(p1, p2, query);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi1 = get(m_vertex_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_vertex_map, *(m_polygon.begin() + (i + 1)));
        A[i] = m_area_2(pi1, pi2, query);
      }

      const auto& pn = get(m_vertex_map, *(m_polygon.begin() + (n - 1)));
      A[n - 1] = m_area_2(pn, p1, query);

      // Initialize weights with areas C following the area notation from [1].
      // Then we multiply them by areas A as in the formula (5) from [1].
      // We also split the loop.
      w[0] = m_area_2(pn, p1, p2);
      for(std::size_t j = 1; j < n - 1; ++j)
        w[0] *= A[j];

      for(std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi0 = get(m_vertex_map, *(m_polygon.begin() + (i - 1)));
        const auto& pi1 = get(m_vertex_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_vertex_map, *(m_polygon.begin() + (i + 1)));
        w[i] = m_area_2(pi0, pi1, pi2);

        for (std::size_t j = 0; j < i - 1; ++j)
          w[i] *= A[j];
        for (std::size_t j = i + 1; j < n; ++j)
          w[i] *= A[j];
      }

      const auto& pm = get(m_vertex_map, *(m_polygon.begin() + (n - 2)));
      w[n - 1] = m_area_2(pm, pn, p1);
      for (std::size_t j = 0; j < n - 2; ++j)
        w[n - 1] *= A[j];

      // Normalize if necessary.
      if (normalize)
        internal::normalize(w);

      // Return weights.
      for (std::size_t i = 0; i < n; ++i)
        *(weights++) = w[i];

      return weights;
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_WACHSPRESS_COORDINATES_2_H
