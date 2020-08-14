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

#ifndef CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_H
#define CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>
#include <CGAL/Weight_interface/internal/polygon_utils.h>

namespace CGAL {
namespace Generalized_weights {

  // [1] Reference: "M. S. Floater, K. Hormann, and G. Kos.
  // A general construction of barycentric coordinates over convex polygons.
  // Advances in Computational Mathematics, 24(1-4):311-331, 2006.".

  // The full weight is computed as
  // \f$w = \frac{C}{A_m A}\f$
  // with notations shown in the figure below. This weight is equal to the
  // `CGAL::Generalized_weights::Authalic_weight`. This weight is a special
  // case of the `CGAL::Generalized_weights::Three_point_family_weight`.
  // \cgalFigureBegin{wachspress_weight, wachspress.svg}
  //   Notation used for the Wachspress weight.
  // \cgalFigureEnd

  /// \cond SKIP_IN_MANUAL
  namespace internal {

    template<typename FT>
    const FT weight(
      const FT A1, const FT A2, const FT C) {

      FT w = FT(0);
      CGAL_assertion(A1 != FT(0) && A2 != FT(0));
      const FT prod = A1 * A2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = C * inv;
      }
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 2D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) wachspress_weight_2(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT C  = internal::area_2(traits, t, r, p);
    return internal::weight(A1, A2, C);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 2D points.

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
  decltype(auto) wachspress_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return wachspress_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 3D points.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2` and `AnalyticTraits_3`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \param traits
    an instance of `GeomTraits`

    \return the computed weight.
  */
  template<typename GeomTraits>
  decltype(auto) wachspress_weight_3(
    const typename GeomTraits::Point_3& q,
    const typename GeomTraits::Point_3& t,
    const typename GeomTraits::Point_3& r,
    const typename GeomTraits::Point_3& p,
    const GeomTraits& traits) {

    using Point_2 = typename GeomTraits::Point_2;
    Point_2 qf, tf, rf, pf;
    internal::flatten(
      traits,
      q,  t,  r,  p,
      qf, tf, rf, pf);
    return wachspress_weight_2(qf, tf, rf, pf, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Wachspress weight for 3D points.

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
  decltype(auto) wachspress_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return wachspress_weight_3(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief 2D Wachspress weights for polygons.

    This class implements 2D Wachspress weights ( \cite cgal:bc:fhk-gcbcocp-06,
    \cite cgal:bc:mlbd-gbcip-02, \cite cgal:bc:w-rfeb-75 ), which can be computed
    at any point inside a strictly convex polygon.

    Wachspress weights are well-defined and non-negative in the closure
    of a strictly convex polygon. The weights are computed analytically.
    See more details in the user manual here.

    \tparam Polygon
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \tparam VertexMap
    must be a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Wachspress_weights_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    using Area_2 = typename GeomTraits::Compute_area_2;
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

      This class implements the behavior of Wachspress weights
      for 2D query points.

      \param polygon
      An instance of `Polygon` with the vertices of a strictly convex polygon.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon`
      to `Point_2`. The default is the identity property map.

      \pre polygon.size() >= 3
      \pre polygon is simple
      \pre polygon is strictly convex
    */
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

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D Wachspress weights.

      This function fills `weights` with 2D Wachspress weights computed at the `query`
      point with respect to the vertices of the input polygon. If `query` belongs to
      the polygon boundary, the returned weights are not defined. You can see the more
      precise version in the package Barycentric Coordinates 2.

      The number of returned weights equals to the number of polygon vertices.

      \tparam OutputIterator
      the dereferenced output iterator type must be convertible to `FT`.

      \param query
      A query point.

      \param w_begin
      The beginning of the destination range with the computed weights.

      \return an output iterator to the element in the destination range,
      one past the last weight stored.
    */
    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator w_begin) {

      const bool normalize = false;
      return optimal_weights(
        query, w_begin, normalize);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<typename OutputIterator>
    OutputIterator operator()(
      const Point_2& query,
      OutputIterator weights,
      const bool normalize) {

      return optimal_weights(
        query, weights, normalize);
    }
    /// \endcond

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

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `CGAL::Generalized_weights::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre polygon.size() >= 3
    \pre polygon is simple
    \pre polygon is strictly convex
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator w_begin,
    const GeomTraits traits) {

    Wachspress_weights_2<PointRange, GeomTraits> wachspress(
      polygon, traits);
    return wachspress(query, w_begin);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `CGAL::Generalized_weights::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre polygon.size() >= 3
    \pre polygon is simple
    \pre polygon is strictly convex
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator w_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_weights_2(
      polygon, query, w_begin, GeomTraits());
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WACHSPRESS_WEIGHTS_H
