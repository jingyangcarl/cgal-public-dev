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

#ifndef CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>
#include <CGAL/Barycentric_coordinates_2/Wachspress_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_coordinates_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `p0 != p1`
  */
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator segment_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits) {

    return internal::linear_coordinates_2(
      p0, p1, query, c_begin, traits);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    must be a point type.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `Kernel_traits<Point_2>::Kernel::FT`.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `p0 != p1`
  */
  template<
  typename Point_2,
  typename OutputIterator>
  OutputIterator segment_coordinates_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& query,
    OutputIterator c_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return segment_coordinates_2(
      p0, p1, query, c_begin, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in a pair.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

    \param query
    A query point.

    \param traits
    An instance of `GeomTraits`.

    \return an `std::pair<GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates.

    \pre `p0 != p1`
  */
  template<typename GeomTraits>
  std::pair<
  typename GeomTraits::FT,
  typename GeomTraits::FT>
  segment_coordinates_in_pair_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits) {

    using FT = typename GeomTraits::FT;
    std::vector<FT> coordinates;
    coordinates.reserve(2);
    internal::linear_coordinates_2(
      p0, p1, query, std::back_inserter(coordinates), traits);
    CGAL_assertion(coordinates.size() == 2);
    return std::make_pair(coordinates[0], coordinates[1]);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in a pair.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    must be a point type.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

    \param query
    A query point.

    \return an `std::pair<FT, FT>` with the computed coordinates, where
    `FT = Kernel_traits<Point_2>::Kernel::FT`.

    \pre `p0 != p1`
  */
  template<typename Point_2>
  std::pair<
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT>
  segment_coordinates_in_pair_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& query) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return segment_coordinates_in_pair_2(
      p0, p1, query, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `area_2(p0, p1, p2) != 0`
  */
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator triangle_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits) {

    return internal::planar_coordinates_2(
      p0, p1, p2, query, c_begin, traits);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    must be a point type.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `Kernel_traits<Point_2>::Kernel::FT`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `area_2(p0, p1, p2) != 0`
  */
  template<
  typename Point_2,
  typename OutputIterator>
  OutputIterator triangle_coordinates_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& p2,
    const Point_2& query,
    OutputIterator c_begin) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return triangle_coordinates_2(
      p0, p1, p2, query, c_begin, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in a tuple.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param traits
    An instance of `GeomTraits`.

    \return an `std::tuple<GeomTraits::FT, GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates.

    \pre `area_2(p0, p1, p2) != 0`
  */
  template<typename GeomTraits>
  std::tuple<
  typename GeomTraits::FT,
  typename GeomTraits::FT,
  typename GeomTraits::FT>
  triangle_coordinates_in_tuple_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits) {

    using FT = typename GeomTraits::FT;
    std::vector<FT> coordinates;
    coordinates.reserve(3);
    internal::planar_coordinates_2(
      p0, p1, p2, query, std::back_inserter(coordinates), traits);
    CGAL_assertion(coordinates.size() == 3);
    return std::make_tuple(coordinates[0], coordinates[1], coordinates[2]);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in a tuple.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    must be a point type.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \return an `std::tuple<FT, FT, FT>` with the computed coordinates, where
    `FT = Kernel_traits<Point_2>::Kernel::FT`.

    \pre `area_2(p0, p1, p2) != 0`
  */
  template<typename Point_2>
  std::tuple<
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT,
  typename Kernel_traits<Point_2>::Kernel::FT>
  triangle_coordinates_in_tuple_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& p2,
    const Point_2& query) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return triangle_coordinates_in_tuple_2(
      p0, p1, p2, query, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D boundary coordinates.

    This function computes boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \tparam VertexMap
    must be a `ReadablePropertyMap` whose key type is `PointRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon`
    to `GeomTraits::Point_2`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polygon boundary.

    \pre `polygon.size() >= 3`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    const auto result =
    internal::locate_wrt_polygon_2(polygon, query, traits, vertex_map);
    auto location = (*result).first;
    auto index = (*result).second;

    if (!result) {
      index = std::size_t(-1);
      location = internal::Query_point_location::UNSPECIFIED;
    }
    return internal::boundary_coordinates_2(
      polygon, query, location, index, c_begin, traits, vertex_map);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D boundary coordinates.

    This function computes 2D boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam Point_2
    must be a point type.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `Kernel_traits<Point_2>::Kernel::FT`.

    \tparam VertexMap
    must be a `ReadablePropertyMap` whose key type is `PointRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon`
    to `Point_2`. The default is the identity property map.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored + the flag indicating whether the
    query point belongs to the polygon boundary.

    \pre `polygon.size() >= 3`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator,
  typename VertexMap = CGAL::Identity_property_map<Point_2> >
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator c_begin,
    const VertexMap vertex_map = VertexMap()) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return boundary_coordinates_2(
      polygon, query, c_begin, GeomTraits(), vertex_map);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator w_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Wachspress_coordinates_2<PointRange, GeomTraits> wachspress(
      polygon, policy, traits);
    return wachspress.weights(query, w_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_coordinates_2` is used.
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

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator w_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_weights_2(
      polygon, query, w_begin, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress coordinates.

    This function computes 2D Wachspress coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator wachspress_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Wachspress_coordinates_2<PointRange, GeomTraits> wachspress(
      polygon, policy, traits);
    return wachspress(query, c_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress coordinates.

    This function computes 2D Wachspress coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_coordinates_2` is used.
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

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator wachspress_coordinates_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator c_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_coordinates_2(
      polygon, query, c_begin, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator w_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Discrete_harmonic_coordinates_2<PointRange, GeomTraits> discrete_harmonic(
      polygon, policy, traits);
    return discrete_harmonic.weights(query, w_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2` is used.
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

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator w_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return discrete_harmonic_weights_2(
      polygon, query, w_begin, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic coordinates.

    This function computes 2D discrete harmonic coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator discrete_harmonic_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Discrete_harmonic_coordinates_2<PointRange, GeomTraits> discrete_harmonic(
      polygon, policy, traits);
    return discrete_harmonic(query, c_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic coordinates.

    This function computes 2D discrete harmonic coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_coordinates_2` is used.
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

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator discrete_harmonic_coordinates_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator c_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return discrete_harmonic_coordinates_2(
      polygon, query, c_begin, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param w_begin
    The beginning of the destination range with the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator w_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Mean_value_coordinates_2<PointRange, GeomTraits> mean_value(
      polygon, policy, traits);
    return mean_value.weights(query, w_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are stored in a destination range
    beginning at 'w_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_coordinates_2` is used.
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

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last weight stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator w_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_weights_2(
      polygon, query, w_begin, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value coordinates.

    This function computes 2D mean value coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_coordinates_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    the dereferenced output iterator type must be convertible to `GeomTraits::FT`.

    \tparam GeomTraits
    must be a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator mean_value_coordinates_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator c_begin,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Mean_value_coordinates_2<PointRange, GeomTraits> mean_value(
      polygon, policy, traits);
    return mean_value(query, c_begin);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value coordinates.

    This function computes 2D mean value coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are stored in a destination range
    beginning at 'c_begin'.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_coordinates_2` is used.
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

    \param c_begin
    The beginning of the destination range with the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy_2`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy_2::DEFAULT`.

    \return an output iterator to the element in the destination range,
    one past the last coordinate stored.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator mean_value_coordinates_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator c_begin,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_coordinates_2(
      polygon, query, c_begin, GeomTraits(), policy);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
