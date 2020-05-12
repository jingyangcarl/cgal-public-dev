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
#include <CGAL/Barycentric_coordinates_2/Wachspress_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/Discrete_harmonic_weights_2.h>
#include <CGAL/Barycentric_coordinates_2/Mean_value_weights_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in `coordinates`.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

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

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in `coordinates`.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    This function infers a traits class from the `Point_2` class.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param p0
    The first vertex of a segment.

    \param p1
    The second vertex of a segment.

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

    \brief computes segment coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the end points `p0` and `p1` of a segment that is one
    coordinate per end point. The coordinates are returned in a pair.

    If `query` does not belong to the line through `p0` and `p1`, it is
    projected onto this line, and only then the coordinates are computed.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

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
    is a point type.

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
    coordinate per point. The coordinates are returned in `coordinates`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

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

    \brief computes triangle coordinates.

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

    \brief computes triangle coordinates.

    This function computes barycentric coordinates at a given `query` point
    with respect to the points `p0`, `p1`, and `p2`, which form a triangle, that is one
    coordinate per point. The coordinates are returned in a tuple.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

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
    is a point type.

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
    coordinate per vertex. The coordinates are returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `PointRange::value_type` and
    value type is `GeomTraits::Point_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

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
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    CGAL_precondition(polygon.size() >= 3);
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

    \brief computes 2D boundary coordinates.

    This function computes 2D boundary barycentric coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    If `query` is at the vertex, the corresponding coordinate is set to one, while
    all other coordinates are zero. If `query` is on the edge, the two corresponding
    coordinates are segment coordinates, while all other coordinates are set to zero.
    If `query` is not on the boundary, all the coordinates are set to zero.

    Internally, `CGAL::Barycentric_coordinates::segment_coordinates_2()` are used.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam Point_2
    is a point type.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `PointRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param vertex_map
    An instance of `VertexMap` that maps a vertex from `polygon`
    to `Point_2`. The default is the identity property map.

    \return an output iterator with the flag indicating whether
    the query point belongs to the polygon boundary.

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
    OutputIterator coordinates,
    const VertexMap vertex_map = VertexMap()) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return boundary_coordinates_2(
      polygon, query, coordinates, GeomTraits(), vertex_map);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Wachspress_weights_2<GeomTraits> wachspress(
      polygon, policy, traits);
    return wachspress(query, weights);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress coordinates.

    This function computes 2D Wachspress coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Wachspress_weights_2<GeomTraits> wachspress(
      polygon, policy, traits);
    return wachspress.coordinates(query, coordinates);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D Wachspress coordinates.

    This function computes 2D Wachspress coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_coordinates_2(
      polygon, query, coordinates, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Discrete_harmonic_weights_2<GeomTraits> discrete_harmonic(
      polygon, policy, traits);
    return discrete_harmonic(query, weights);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return discrete_harmonic_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic coordinates.

    This function computes 2D discrete harmonic coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Discrete_harmonic_weights_2<GeomTraits> discrete_harmonic(
      polygon, policy, traits);
    return discrete_harmonic.coordinates(query, coordinates);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D discrete harmonic coordinates.

    This function computes 2D discrete harmonic coordinates at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return discrete_harmonic_coordinates_2(
      polygon, query, coordinates, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Mean_value_weights_2<GeomTraits> mean_value(
      polygon, policy, traits);
    return mean_value(query, weights);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator weights,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value coordinates.

    This function computes 2D mean value coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const GeomTraits traits,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    Mean_value_weights_2<GeomTraits> mean_value(
      polygon, policy, traits);
    return mean_value.coordinates(query, coordinates);
  }

  /*!
    \ingroup PkgBarycentricCoordinates2RefFunctions

    \brief computes 2D mean value coordinates.

    This function computes 2D mean value coordinates at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    coordinate per vertex. The coordinates are returned in `coordinates`.

    Internally, the class `CGAL::Barycentric_coordinates::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param coordinates
    An output iterator that stores the computed coordinates.

    \param policy
    One of the `CGAL::Barycentric_coordinates::Computation_policy`.
    The default is `CGAL::Barycentric_coordinates::Computation_policy::DEFAULT`.

    \return an output iterator.

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
    OutputIterator coordinates,
    const Computation_policy policy =
    Computation_policy::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_coordinates_2(
      polygon, query, coordinates, GeomTraits(), policy);
  }

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ANALYTIC_COORDINATES_2_H
