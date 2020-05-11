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

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in `coordinates`.

    If the `query` point does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

    \param p0
    The source vertex of a segment.

    \param p1
    The target vertex of a segment.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an output iterator.

    \pre `p0 != p1`
  */
  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator segment_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits) {

    return internal::linear_coordinates_2(
      p0, p1, query, coordinates, traits);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in `coordinates`.

    If the `query` point does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param p0
    The source vertex of a segment.

    \param p1
    The target vertex of a segment.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \return an output iterator.

    \pre `p0 != p1`
  */
  template<
  typename Point_2,
  typename OutputIterator>
  OutputIterator segment_coordinates_2(
    const Point_2& p0,
    const Point_2& p1,
    const Point_2& query,
    OutputIterator coordinates) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return segment_coordinates_2(
      p0, p1, query, coordinates, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point.

    If the `query` point does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

    \param p0
    The source vertex of a segment.

    \param p1
    The target vertex of a segment.

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
    return std::make_pair(coordinates[0], coordinates[1]);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point.

    If the `query` point does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    is a point type.

    \param p0
    The source vertex of a segment.

    \param p1
    The target vertex of a segment.

    \param query
    A query point.

    \return an `std::pair<GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates.

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

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in `coordinates`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const GeomTraits traits) {

    return internal::planar_coordinates_2(
      p0, p1, p2, query, coordinates, traits);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in `coordinates`.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \return an output iterator.

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
    OutputIterator coordinates) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return triangle_coordinates_2(
      p0, p1, p2, query, coordinates, GeomTraits());
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

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
    return std::make_tuple(coordinates[0], coordinates[1], coordinates[2]);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    is a point type.

    \param p0
    The first vertex of a triangle.

    \param p1
    The second vertex of a triangle.

    \param p2
    The third vertex of a triangle.

    \param query
    A query point.

    \return an `std::tuple<GeomTraits::FT, GeomTraits::FT, GeomTraits::FT>`
    with the computed coordinates.

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

    This function takes a `query` point and computes boundary barycentric
    coordinates at this point with respect to the vertices of a given `polygon`.
    These coordinates are then returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `GeomTraits::Point_2`.

    \param polygon
    An instance of `Polygon` with the vertices of a simple polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon`
    to `GeomTraits::Point_2`.

    \return an output iterator with the flag indicating whether
    the query point belongs to the polygon boundary.
  */
  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using Point_2 = typename GeomTraits::Point_2;
    std::vector<Point_2> poly;
    poly.reserve(polygon.size());
    for (const auto& item : polygon)
      poly.push_back(get(vertex_map, item));

    const auto result =
    internal::locate_wrt_polygon_2(poly, query, traits);
    auto location = (*result).first;
    auto index    = (*result).second;

    if (!result) {
      location = internal::Query_point_location::UNSPECIFIED;
      index = std::size_t(-1);
    }

    return internal::boundary_coordinates_2(
      poly, query, location, index, coordinates, traits);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function takes a `query` point and computes boundary barycentric
    coordinates at this point with respect to the vertices of a given `polygon`.
    These coordinates are then returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    This function infers a traits class from the `Point_2` class.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `Polygon` with the vertices of a simple polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon`
    to `Point_2`. The default is the identity property map.

    \return an output iterator with the flag indicating whether
    the query point belongs to the polygon boundary.
  */
  template<
  typename Polygon,
  typename Point_2,
  typename OutputIterator,
  typename VertexMap = CGAL::Identity_property_map<
  typename Kernel_traits<Point_2>::Kernel::Point_2> >
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const Polygon& polygon,
    const Point_2& query,
    OutputIterator coordinates,
    const VertexMap vertex_map = VertexMap()) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return boundary_coordinates_2(
      polygon, query, coordinates, GeomTraits(), vertex_map);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function takes a `query` point and computes the chosen barycentric
    `weights` at this point with respect to the given `vertices`. These weights
    are then normalized and returned in `coordinates`.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam Weights
    is a model of `CGAL::Barycentric_coordinates::AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `CGAL::Barycentric_coordinates::BarycentricTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param vertices
    An instance of `VertexRange` with vertices, e.g. polygon vertices.

    \param query
    A query point.

    \param weights
    An instance of `Weights` that computes the corresponding
    barycentric weights at the given query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `vertices`
    to `GeomTraits::Point_2`.

    \return an output iterator.
  */
  template<
  typename VertexRange,
  typename Point_2,
  typename Weights,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  OutputIterator analytic_coordinates_2(
    const VertexRange& vertices,
    const Point_2& query,
    Weights& weights,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using FT = typename GeomTraits::FT;

    std::vector<FT> b;
    b.reserve(vertices.size());

    weights(
      vertices, query, std::back_inserter(b), traits, vertex_map);
    internal::normalize(b);
    for (const auto& value : b)
      *(coordinates++) = value;
    return coordinates;
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    This function takes a `query` point and computes the chosen barycentric
    `weights` at this point with respect to the given `vertices`. These weights
    are then normalized and returned in `coordinates`.

    This function infers a traits class from the `Point_2` class.

    \tparam VertexRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam Weights
    is a model of `CGAL::Barycentric_coordinates::AnalyticWeights_2`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `GeomTraits::Point_2`. The default is `CGAL::Identity_property_map`.

    \param vertices
    An instance of `VertexRange` with vertices, e.g. polygon vertices.

    \param query
    A query point.

    \param weights
    An instance of `Weights` that computes the corresponding
    barycentric weights at the given query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `vertices`
    to `GeomTraits::Point_2`. The default is the identity property map.

    \return an output iterator.
  */
  template<
  typename VertexRange,
  typename Point_2,
  typename Weights,
  typename OutputIterator,
  typename VertexMap = CGAL::Identity_property_map<typename Kernel_traits<Point_2>::Kernel::Point_2> >
  OutputIterator analytic_coordinates_2(
    const VertexRange& vertices,
    const Point_2& query,
    Weights& weights,
    OutputIterator coordinates,
    const VertexMap vertex_map = VertexMap()) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return analytic_coordinates_2(
      vertices, query, weights, coordinates, GeomTraits(), vertex_map);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
