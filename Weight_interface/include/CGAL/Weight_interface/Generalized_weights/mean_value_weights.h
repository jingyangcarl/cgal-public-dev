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

#ifndef CGAL_GENERALIZED_MEAN_VALUE_WEIGHTS_H
#define CGAL_GENERALIZED_MEAN_VALUE_WEIGHTS_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>
#include <CGAL/Weight_interface/internal/polygon_utils.h>

namespace CGAL {
namespace Generalized_weights {

  // [1] Reference: "K. Hormann and M. Floater.
  // Mean value coordinates for arbitrary planar polygons.
  // ACM Transactions on Graphics, 25(4):1424-1441, 2006.".

  // The full weight is computed as
  // \f$w = \pm 2 \sqrt{\frac{2 (r_m r_p - D)}{(r r_m + D_m)(r r_p + D_p)}}\f$,
  // with notations shown in the figure below and dot products
  // \f$D_m = (v_j - q) \cdot (v_m - q)\f$,
  // \f$D_p = (v_j - q) \cdot (v_p - q)\f$, and
  // \f$D   = (v_m - q) \cdot (v_p - q)\f$.
  // The \f$\pm\f$ sign is a sign of the weight that depends on the configuration.
  // This weight is equal to the `CGAL::Generalized_weights::Tangent_weight`.
  // This weight is a special case of the `CGAL::Generalized_weights::Three_point_family_weight`.
  // \cgalFigureBegin{mean_value_weight, mean_value.svg}
  //   Notation used for the mean value weight.
  // \cgalFigureEnd

  /// \cond SKIP_IN_MANUAL
  namespace mean_value_ns {

    template<typename FT>
    const FT sign_of_weight(
      const FT A1, const FT A2, const FT B) {

      if (A1 > FT(0) && A2 > FT(0) && B <= FT(0)) return  FT(1);
      if (A1 < FT(0) && A2 < FT(0) && B >= FT(0)) return -FT(1);
      if (B  > FT(0)) return  FT(1);
      if (B  < FT(0)) return -FT(1);
      return FT(0);
    }

    template<typename GeomTraits>
    const typename GeomTraits::FT weight(
      const GeomTraits& traits,
      const typename GeomTraits::FT r1,
      const typename GeomTraits::FT r2,
      const typename GeomTraits::FT r3,
      const typename GeomTraits::FT D1,
      const typename GeomTraits::FT D2,
      const typename GeomTraits::FT D,
      const typename GeomTraits::FT sign) {

      using FT = typename GeomTraits::FT;
      using Get_sqrt = internal::Get_sqrt<GeomTraits>;
      const auto sqrt = Get_sqrt::sqrt_object(traits);

      const FT P1 = r1 * r2 + D1;
      const FT P2 = r2 * r3 + D2;

      FT w = FT(0);
      CGAL_assertion(P1 != FT(0) && P2 != FT(0));
      const FT prod = P1 * P2;
      if (prod != FT(0)) {
        const FT inv = FT(1) / prod;
        w = FT(2) * (r1 * r3 - D) * inv;
        CGAL_assertion(w >= FT(0));
        w = sqrt(w);
      }
      w *= FT(2); w *= sign;
      return w;
    }
  }
  /// \endcond

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the mean value weight for 2D points.

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
  decltype(auto) mean_value_weight_2(
    const typename GeomTraits::Point_2& q,
    const typename GeomTraits::Point_2& t,
    const typename GeomTraits::Point_2& r,
    const typename GeomTraits::Point_2& p,
    const GeomTraits& traits) {

    using FT = typename GeomTraits::FT;
    const auto dot_product_2 =
      traits.compute_scalar_product_2_object();
    const auto construct_vector_2 =
      traits.construct_vector_2_object();

    const auto v1 = construct_vector_2(q, t);
    const auto v2 = construct_vector_2(q, r);
    const auto v3 = construct_vector_2(q, p);

    const FT l1 = internal::length_2(traits, v1);
    const FT l2 = internal::length_2(traits, v2);
    const FT l3 = internal::length_2(traits, v3);

    const FT D1 = dot_product_2(v1, v2);
    const FT D2 = dot_product_2(v2, v3);
    const FT D  = dot_product_2(v1, v3);

    const FT A1 = internal::area_2(traits, r, q, t);
    const FT A2 = internal::area_2(traits, p, q, r);
    const FT B  = internal::area_2(traits, p, q, t);

    const FT sign = mean_value_ns::sign_of_weight(A1, A2, B);
    return mean_value_ns::weight(
      traits, l1, l2, l3, D1, D2, D, sign);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPoints

    \brief computes the mean value weight for 2D points.

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
  decltype(auto) mean_value_weight_2(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    const GeomTraits traits;
    return mean_value_weight_2(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the mean value weight for 3D points.

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
  decltype(auto) mean_value_weight_3(
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
    return mean_value_weight_2(qf, tf, rf, pf, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights3DPoints

    \brief computes the mean value weight for 3D points.

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
  decltype(auto) mean_value_weight_3(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) {

    using GeomTraits = typename Kernel_traits<Point_3>::Kernel;
    const GeomTraits traits;
    return mean_value_weight_3(q, t, r, p, traits);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPolygons

    \brief 2D mean value weights for polygons.

    This class implements 2D mean value weights ( \cite cgal:bc:hf-mvcapp-06,
    \cite cgal:bc:fhk-gcbcocp-06, \cite cgal:f-mvc-03 ), which can be computed
    at any point in the plane.

    Mean value weights are well-defined everywhere in the plane and are
    non-negative in the kernel of a star-shaped polygon. The weights are
    computed analytically. See more details in the user manual here.

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
  class Mean_value_weights_2 {

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

      \pre polygon.size() >= 3
      \pre polygon is simple
    */
    Mean_value_weights_2(
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
      \brief computes 2D mean value weights.

      This function fills `weights` with 2D mean value weights computed at the `query`
      point with respect to the vertices of the input polygon. If `query` belongs to
      the polygon boundary, the weights are not well-defined. You can see the more
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
      return operator()(query, w_begin, normalize);
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
      // in the Figure 10 from [1]. Split the loop to make this computation faster.
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
        CGAL_assertion((r[i] * r[i + 1] + D[i]) != FT(0));
        t[i] = FT(2) * A[i] / (r[i] * r[i + 1] + D[i]);
      }

      CGAL_assertion((r[n - 1] * r[0] + D[n - 1]) != FT(0));
      t[n - 1] = FT(2) * A[n - 1] / (r[n - 1] * r[0] + D[n - 1]);

      // Compute mean value weights using the same pseudo-code as before.
      CGAL_assertion(r[0] != FT(0));
      w[0] = FT(2) * (t[n - 1] + t[0]) / r[0];

      for (std::size_t i = 1; i < n - 1; ++i) {
        CGAL_assertion(r[i] != FT(0));
        w[i] = FT(2) * (t[i - 1] + t[i]) / r[i];
      }

      CGAL_assertion(r[n - 1] != FT(0));
      w[n - 1] = FT(2) * (t[n - 2] + t[n - 1]) / r[n - 1];

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
    \ingroup PkgWeightInterfaceRefWeights2DPolygons

    \brief computes 2D mean value weights for polygons.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `CGAL::Generalized_weights::Mean_value_weights_2` is used.
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
    An instance of `PointRange` with 2D points, which form a simple polygon.

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
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator w_begin,
    const GeomTraits traits) {

    Mean_value_weights_2<PointRange, GeomTraits> mean_value(
      polygon, traits);
    return mean_value(query, w_begin);
  }

  /*!
    \ingroup PkgWeightInterfaceRefWeights2DPolygons

    \brief computes 2D mean value weights for polygons.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at `w_begin`.

    Internally, the class `CGAL::Generalized_weights::Mean_value_weights_2` is used.
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
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre polygon.size() >= 3
    \pre polygon is simple
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator w_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_weights_2(
      polygon, query, w_begin, GeomTraits());
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_MEAN_VALUE_WEIGHTS_H
