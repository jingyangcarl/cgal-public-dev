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

#ifndef CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H
#define CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/Generalized_weights_2/internal/utils_2.h>

// [1] Reference: "K. Hormann and M. Floater.
// Mean value coordinates for arbitrary planar polygons.
// ACM Transactions on Graphics, 25(4):1424-1441, 2006.".

// [2] Reference: "M. S. Floater.
// Wachspress and mean value coordinates.
// Proceedings of the 14th International Conference on Approximation Theory.
// G. Fasshauer and L. L. Schumaker (eds.)."

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2D

    \brief 2D mean value weight.

    This class implements 2D mean value weight ( \cite cgal:bc:hf-mvcapp-06,
    \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ).

    \tparam Polygon
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits
    is a model of `AnalyticTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \cgalModels `AnalyticWeight_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Mean_value_weight_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    using Vector_2 = typename GeomTraits::Vector_2;
    using Area_2 = typename GeomTraits::Compute_area_2;
    using Squared_length_2 = typename GeomTraits::Compute_squared_length_2;
    using Scalar_product_2 = typename GeomTraits::Compute_scalar_product_2;
    using Get_sqrt = internal::Get_sqrt<GeomTraits>;
    using Sqrt = typename Get_sqrt::Sqrt;
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

      This class implements the behavior of mean value weights
      for 2D query points.

      \param polygon
      An instance of `Polygon` with the vertices of a simple polygon.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon`
      to `Point_2`. The default is the identity property map.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
    */
    Mean_value_weight_2(
      const Polygon& polygon,
      const GeomTraits traits = GeomTraits(),
      const VertexMap vertex_map = VertexMap()) :
    m_polygon(polygon),
    m_traits(traits),
    m_vertex_map(vertex_map),
    m_area_2(m_traits.compute_area_2_object()),
    m_squared_length_2(m_traits.compute_squared_length_2_object()),
    m_scalar_product_2(m_traits.compute_scalar_product_2_object()),
    m_sqrt(Get_sqrt::sqrt_object(m_traits))  {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, vertex_map));
      resize();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D mean value weight.

      This function fills `weights` with 2D mean value weights computed at the `query`
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

  private:

    // Fields.
    const Polygon& m_polygon;
    const Computation_policy_2 m_computation_policy;
    const GeomTraits m_traits;
    const VertexMap m_vertex_map;
    const Area_2 m_area_2;
    const Squared_length_2 m_squared_length_2;
    const Scalar_product_2 m_scalar_product_2;
    const Sqrt m_sqrt;

    std::vector<Vector_2> s;
    std::vector<FT> r;
    std::vector<FT> A;
    std::vector<FT> D;
    std::vector<FT> t;
    std::vector<FT> w;

    // Functions.
    void resize() {
      s.resize(m_polygon.size());
      r.resize(m_polygon.size());
      A.resize(m_polygon.size());
      D.resize(m_polygon.size());
      t.resize(m_polygon.size());
      w.resize(m_polygon.size());
    }

    template<typename OutputIterator>
    OutputIterator optimal_weights(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      // Get the number of vertices in the polygon.
      const std::size_t n = m_polygon.size();

      // Compute vectors s following the pseudo-code in the Figure 10 from [1].
      for (std::size_t i = 0; i < n; ++i) {
        const auto& pi = get(m_vertex_map, *(m_polygon.begin() + i));
        s[i] = pi - query;
      }

      // Compute lengths r, areas A, and dot products D following the pseudo-code
      // in the Figure 10 from [1].
      // Split the loop to make this computation faster.
      const auto& p1 = get(m_vertex_map, *(m_polygon.begin() + 0));
      const auto& p2 = get(m_vertex_map, *(m_polygon.begin() + 1));

      r[0] = m_sqrt(m_squared_length_2(s[0]));
      A[0] = m_area_2(p1, p2, query);
      D[0] = m_scalar_product_2(s[0], s[1]);

      for (std::size_t i = 1; i < n - 1; ++i) {
        const auto& pi1 = get(m_vertex_map, *(m_polygon.begin() + (i + 0)));
        const auto& pi2 = get(m_vertex_map, *(m_polygon.begin() + (i + 1)));

        r[i] = m_sqrt(m_squared_length_2(s[i]));
        A[i] = m_area_2(pi1, pi2, query);
        D[i] = m_scalar_product_2(s[i], s[i + 1]);
      }

      const auto& pn = get(m_vertex_map, *(m_polygon.begin() + (n - 1)));
      r[n - 1] = m_sqrt(m_squared_length_2(s[n - 1]));
      A[n - 1] = m_area_2(pn, p1, query);
      D[n - 1] = m_scalar_product_2(s[n - 1], s[0]);

      // Compute intermediate values t using the formulas from slide 19 here
      // - http://www.inf.usi.ch/hormann/nsfworkshop/presentations/Hormann.pdf
      for (std::size_t i = 0; i < n - 1; ++i) {
        CGAL_precondition((r[i] * r[i + 1] + D[i]) != FT(0));
        t[i] = A[i] / (r[i] * r[i + 1] + D[i]);
      }

      CGAL_precondition((r[n - 1] * r[0] + D[n - 1]) != FT(0));
      t[n - 1] = A[n - 1] / (r[n - 1] * r[0] + D[n - 1]);

      // Compute mean value weights using the same pseudo-code as before.
      CGAL_precondition(r[0] != FT(0));
      w[0] = (t[n - 1] + t[0]) / r[0];

      for (std::size_t i = 1; i < n - 1; ++i) {
        CGAL_precondition(r[i] != FT(0));
        w[i] = (t[i - 1] + t[i]) / r[i];
      }

      CGAL_precondition(r[n - 1] != FT(0));
      w[n - 1] = (t[n - 2] + t[n - 1]) / r[n - 1];

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

#endif // CGAL_GENERALIZED_MEAN_VALUE_WEIGHT_2_H
