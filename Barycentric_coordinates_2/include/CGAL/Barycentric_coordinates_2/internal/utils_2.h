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

#ifndef CGAL_BARYCENTRIC_UTILS_2_H
#define CGAL_BARYCENTRIC_UTILS_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// STL includes.
#include <set>
#include <map>
#include <list>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <tuple>

// Boost headers.
#include <boost/mpl/has_xxx.hpp>
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/barycenter.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

namespace CGAL {
namespace Barycentric_coordinates {
namespace internal {

  enum class Edge_case {

    EXTERIOR = 0, // exterior part of the polygon
    BOUNDARY = 1, // boundary part of the polygon
    INTERIOR = 2  // interior part of the polygon
  };

  enum class Query_point_location {

    // Query point is located at the vertex of the polygon.
    ON_VERTEX = 0,

    // Query point is located on the edge of the polygon.
    ON_EDGE = 1,

    // Query point is located in the polygon's interior.
    ON_BOUNDED_SIDE = 2,

    // Query point is located in the polygon's exterior.
    ON_UNBOUNDED_SIDE = 3,

    // Location is unspecified. Leads to all coordinates being set to zero.
    UNSPECIFIED = 4
  };

  enum class Polygon_type {

    // Concave polygon = non-convex polygon.
    CONCAVE = 0,

    // This is a convex polygon with collinear vertices.
    WEAKLY_CONVEX = 1,

    // This is a convex polygon without collinear vertices.
    STRICTLY_CONVEX = 2
  };

  template<typename GeomTraits>
  class Default_sqrt {

  private:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;

  public:
    FT operator()(const FT value) const {
      return static_cast<FT>(
        CGAL::sqrt(CGAL::to_double(CGAL::abs(value))));
    }
  };

  BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Sqrt, Sqrt, false)

  // Case: do_not_use_default = false.
  template<typename GeomTraits,
  bool do_not_use_default = Has_nested_type_Sqrt<GeomTraits>::value>
  class Get_sqrt {

  public:
    using Traits = GeomTraits;
    using Sqrt = Default_sqrt<Traits>;

    static Sqrt sqrt_object(const Traits& ) {
      return Sqrt();
    }
  };

  // Case: do_not_use_default = true.
  template<typename GeomTraits>
  class Get_sqrt<GeomTraits, true> {

  public:
    using Traits = GeomTraits;
    using Sqrt = typename Traits::Sqrt;

    static Sqrt sqrt_object(const Traits& traits) {
      return traits.sqrt_object();
    }
  };

  template<typename OutputIterator>
  void get_default(
    const std::size_t n,
    OutputIterator output) {

    for (std::size_t i = 0; i < n; ++i)
      *(output++) = 0;
  }

  template<typename FT>
  void normalize(std::vector<FT>& values) {

    FT sum = FT(0);
    for (const FT value : values)
      sum += value;

    CGAL_assertion(sum != FT(0));
    if (sum == FT(0))
      return;

    const FT inv_sum = FT(1) / sum;
    for (FT& value : values)
      value *= inv_sum;
  }

  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator linear_coordinates_2(
    const typename GeomTraits::Point_2& source,
    const typename GeomTraits::Point_2& target,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits) {

    CGAL_precondition(source != target);
    if (source == target) {
      get_default(2, coordinates);
      return coordinates;
    }

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto scalar_product_2   = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();

    // Project point onto segment.
    const FT opposite_scalar_product =
    scalar_product_2(query - target, source - target);

    // Compute coordinates.
    CGAL_assertion(source != target);
    const FT b0 = opposite_scalar_product / squared_distance_2(source, target);
    const FT b1 = FT(1) - b0;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;

    return coordinates;
  }

  template<
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator planar_coordinates_2(
    const typename GeomTraits::Point_2& p0,
    const typename GeomTraits::Point_2& p1,
    const typename GeomTraits::Point_2& p2,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits) {

    // Number type.
    using FT = typename GeomTraits::FT;

    // Functions.
    const auto area_2 = traits.compute_area_2_object();
    const FT total_area = area_2(p0, p1, p2);

    CGAL_precondition(total_area != FT(0));
    if (total_area == FT(0)) {
      get_default(3, coordinates);
      return coordinates;
    }

    // Compute some related sub-areas.
    const FT A1 = area_2(p1, p2, query);
    const FT A2 = area_2(p2, p0, query);

    // Compute the inverted total area of the triangle.
    CGAL_assertion(total_area != FT(0));
    const FT inverted_total_area = FT(1) / total_area;

    // Compute coordinates.
    const FT b0 = A1 * inverted_total_area;
    const FT b1 = A2 * inverted_total_area;
    const FT b2 = FT(1) - b0 - b1;

    // Return coordinates.
    *(coordinates++) = b0;
    *(coordinates++) = b1;
    *(coordinates++) = b2;

    return coordinates;
  }

  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  boost::optional< std::pair<Query_point_location, std::size_t> >
  get_edge_index_approximate(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using FT = typename GeomTraits::FT;
    using Vector_2 = typename GeomTraits::Vector_2;

    const auto cross_product_2 = traits.compute_determinant_2_object();
    const auto scalar_product_2 = traits.compute_scalar_product_2_object();
    const auto squared_distance_2 = traits.compute_squared_distance_2_object();

    CGAL_precondition(polygon.size() >= 3);
    const std::size_t n = polygon.size();

    const FT half = FT(1) / FT(2);
    const FT tolerance = FT(1) / FT(100000);
    const FT sq_tolerance = tolerance * tolerance;

    for (std::size_t i = 0; i < n; ++i) {
      const auto& p1 = get(vertex_map, *(polygon.begin() + i));

      const FT sq_r = squared_distance_2(query, p1);
      if (sq_r < sq_tolerance)
        return std::make_pair(Query_point_location::ON_VERTEX, i);

      const std::size_t ip = (i + 1) % n;
      const auto& p2 = get(vertex_map, *(polygon.begin() + ip));

      const Vector_2 s1 = Vector_2(query, p1);
      const Vector_2 s2 = Vector_2(query, p2);

      const FT A = half * cross_product_2(s1, s2);
      const FT D = scalar_product_2(s1, s2);

      if (CGAL::abs(A) < tolerance && D < FT(0))
        return std::make_pair(Query_point_location::ON_EDGE, i);
    }
    return boost::none;
  }

  // Why this one does not work for harmonic coordinates?
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  boost::optional< std::pair<Query_point_location, std::size_t> >
  get_edge_index_exact(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    const auto collinear_2 = traits.collinear_2_object();
    const auto collinear_are_ordered_along_line_2 =
      traits.collinear_are_ordered_along_line_2_object();
    CGAL_precondition(polygon.size() >= 3);

    const std::size_t n = polygon.size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto& p1 = get(vertex_map, *(polygon.begin() + i));

      if (p1 == query)
        return std::make_pair(Query_point_location::ON_VERTEX, i);

      const std::size_t ip = (i + 1) % n;
      const auto& p2 = get(vertex_map, *(polygon.begin() + ip));

      if (
        collinear_2(p1, p2, query) &&
        collinear_are_ordered_along_line_2(p1, query, p2)) {

        return std::make_pair(Query_point_location::ON_EDGE, i);
      }
    }
    return boost::none;
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  class Point_2,
  class Orientation_2,
  class CompareX_2>
  int which_side_in_slab(
    const Point_2& query,
    const Point_2& low,
    const Point_2& high,
    const Orientation_2& orientation_2,
    const CompareX_2& compare_x_2) {

    const auto low_x_comp_res = compare_x_2(query, low);
    const auto high_x_comp_res = compare_x_2(query, high);
    if (low_x_comp_res == CGAL::SMALLER) {
      if (high_x_comp_res == CGAL::SMALLER)
          return -1;
    } else {
      switch (high_x_comp_res) {
        case CGAL::LARGER: return 1;
        case CGAL::SMALLER: break;
        case CGAL::EQUAL: return (low_x_comp_res == CGAL::EQUAL) ? 0 : 1;
      }
    }
    switch (orientation_2(low, query, high)) {
      case CGAL::LEFT_TURN:  return  1;
      case CGAL::RIGHT_TURN: return -1;
      default: return 0;
    }
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  Edge_case bounded_side_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    const auto first = polygon.begin();
    const auto last  = polygon.end();

    auto curr = first;
    if (curr == last)
      return Edge_case::EXTERIOR;

    auto next = curr; ++next;
    if (next == last)
      return Edge_case::EXTERIOR;

    const auto compare_x_2 = traits.compare_x_2_object();
    const auto compare_y_2 = traits.compare_y_2_object();
    const auto orientation_2 = traits.orientation_2_object();

    bool is_inside = false;
    auto curr_y_comp_res = compare_y_2(get(vertex_map, *curr), query);

    // Check if the segment (curr, next) intersects
    // the ray { (t, query.y()) | t >= query.x() }.
    do {
      const auto& currp = get(vertex_map, *curr);
      const auto& nextp = get(vertex_map, *next);

      auto next_y_comp_res = compare_y_2(nextp, query);
      switch (curr_y_comp_res) {
        case CGAL::SMALLER:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              break;
            case CGAL::EQUAL:
              switch (compare_x_2(query, nextp)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case CGAL::EQUAL:   return Edge_case::BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::LARGER:
              switch (which_side_in_slab(
                query, currp, nextp, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return Edge_case::BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::EQUAL:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              switch (compare_x_2(query, currp)) {
                case CGAL::SMALLER: is_inside = !is_inside; break;
                case CGAL::EQUAL:   return Edge_case::BOUNDARY;
                case CGAL::LARGER:  break;
              }
              break;
            case CGAL::EQUAL:
              switch (compare_x_2(query, currp)) {
                case CGAL::SMALLER:
                  if (compare_x_2(query, nextp) != CGAL::SMALLER)
                      return Edge_case::BOUNDARY;
                  break;
                case CGAL::EQUAL: return Edge_case::BOUNDARY;
                case CGAL::LARGER:
                  if (compare_x_2(query, nextp) != CGAL::LARGER)
                      return Edge_case::BOUNDARY;
                  break;
              }
              break;
            case CGAL::LARGER:
              if (compare_x_2(query, currp) == CGAL::EQUAL) {
                return Edge_case::BOUNDARY;
              }
              break;
          }
          break;
        case CGAL::LARGER:
          switch (next_y_comp_res) {
            case CGAL::SMALLER:
              switch (which_side_in_slab(
                query, nextp, currp, orientation_2, compare_x_2)) {
                case -1: is_inside = !is_inside; break;
                case  0: return Edge_case::BOUNDARY;
              }
              break;
            case CGAL::EQUAL:
              if (compare_x_2(query, nextp) == CGAL::EQUAL) {
                return Edge_case::BOUNDARY;
              }
              break;
            case CGAL::LARGER:
              break;
          }
          break;
      }

      curr = next;
      curr_y_comp_res = next_y_comp_res;
      ++next;
      if (next == last) next = first;
    }
    while (curr != first);
    return is_inside ? Edge_case::INTERIOR : Edge_case::EXTERIOR;
  }

  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  boost::optional< std::pair<Query_point_location, std::size_t> >
  locate_wrt_polygon_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    const Edge_case type = bounded_side_2(
      polygon, query, traits, vertex_map);

    // Locate point with respect to different polygon locations.
    switch (type) {
      case Edge_case::INTERIOR:
        return std::make_pair(Query_point_location::ON_BOUNDED_SIDE, std::size_t(-1));
      case Edge_case::EXTERIOR:
        return std::make_pair(Query_point_location::ON_UNBOUNDED_SIDE, std::size_t(-1));
      case Edge_case::BOUNDARY:
        return get_edge_index_exact(polygon, query, traits, vertex_map);
      default:
        return std::make_pair(Query_point_location::UNSPECIFIED, std::size_t(-1));
    }
    return boost::none;
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  bool is_convex_2(
    const Polygon& polygon,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    auto first = polygon.begin();
    const auto last  = polygon.end();

    auto prev = first;
    if (prev == last) return true;

    auto curr = prev; ++curr;
    if (curr == last) return true;

    auto next = curr; ++next;
    if (next == last) return true;

    const auto equal_2 = traits.equal_2_object();

    while (equal_2(
      get(vertex_map, *prev),
      get(vertex_map, *curr))) {
      curr = next; ++next;
      if (next == last) return true;
    }

    const auto less_xy_2 = traits.less_xy_2_object();
    const auto orientation_2 = traits.orientation_2_object();

    bool has_clockwise_triplets = false;
    bool has_counterclockwise_triplets = false;
    bool order = less_xy_2(
      get(vertex_map, *prev), get(vertex_map, *curr));
    int num_order_changes = 0;

    do {
    switch_orient:
      switch (orientation_2(
        get(vertex_map, *prev),
        get(vertex_map, *curr),
        get(vertex_map, *next))) {

        case CGAL::CLOCKWISE:
          has_clockwise_triplets = true;
          break;
        case CGAL::COUNTERCLOCKWISE:
          has_counterclockwise_triplets = true;
          break;
        case CGAL::ZERO: {
          if (equal_2(
            get(vertex_map, *curr),
            get(vertex_map, *next))) {

            if (next == first) {
              first = curr;
            }
            ++next;
            if (next == last)
              next = first;
            goto switch_orient;
          }
          break;
        }
      }

      const bool new_order = less_xy_2(
        get(vertex_map, *curr),
        get(vertex_map, *next));
      if (order != new_order) num_order_changes++;
      if (num_order_changes > 2)
        return false;
      if (has_clockwise_triplets && has_counterclockwise_triplets)
        return false;

      prev = curr;
      curr = next;
      ++next;
      if (next == last) next = first;
      order = new_order;
    } while (prev != first);
    return true;
  }

  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  Polygon_type polygon_type_2(
    const Polygon& polygon,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using Point_2 = typename GeomTraits::Point_2;
    const auto collinear_2 = traits.collinear_2_object();
    CGAL_precondition(polygon.size() >= 3);

    // First, test the polygon on convexity.
    if (is_convex_2(polygon, traits, vertex_map)) {

      // Test all the consequent triplets of the polygon vertices on collinearity.
      // In case we find at least one, return WEAKLY_CONVEX polygon.
      const std::size_t n = polygon.size();
      for (std::size_t i = 0; i < n; ++i) {
        const auto& p1 = get(vertex_map, *(polygon.begin() + i));

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        const auto& p0 = get(vertex_map, *(polygon.begin() + im));
        const auto& p2 = get(vertex_map, *(polygon.begin() + ip));

        if (collinear_2(p0, p1, p2))
          return Polygon_type::WEAKLY_CONVEX;
      }
      // Otherwise, return STRICTLY_CONVEX polygon.
      return Polygon_type::STRICTLY_CONVEX;
    }
    // Otherwise, return CONCAVE polygon.
    return Polygon_type::CONCAVE;
  }

  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  std::pair<OutputIterator, bool> coordinates_on_last_edge_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using FT = typename GeomTraits::FT;
    const std::size_t n = polygon.size();

    std::vector<FT> b;
    b.reserve(2);

    const std::size_t isource = n - 1;
    const std::size_t itarget = 0;

    const auto& source = get(vertex_map, *(polygon.begin() + isource));
    const auto& target = get(vertex_map, *(polygon.begin() + itarget));

    linear_coordinates_2(
      source, target, query, std::back_inserter(b), traits);
    *(coordinates++) = b[1];
    for (std::size_t i = 1; i < n - 1; ++i)
      *(coordinates++) = FT(0);
    *(coordinates++) = b[0];

    return std::make_pair(coordinates, true);
  }

  template<
  typename Polygon,
  typename OutputIterator,
  typename GeomTraits,
  typename VertexMap>
  std::pair<OutputIterator, bool> boundary_coordinates_2(
    const Polygon& polygon,
    const typename GeomTraits::Point_2& query,
    const Query_point_location location,
    const std::size_t index,
    OutputIterator coordinates,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    using FT = typename GeomTraits::FT;
    const std::size_t n = polygon.size();

    // Compute coordinates with respect to the query point location.
    switch (location) {

      case Query_point_location::ON_VERTEX: {
        CGAL_assertion(index >= 0 && index < n);

        for (std::size_t i = 0; i < n; ++i)
          if (i == index)
            *(coordinates++) = FT(1);
          else
            *(coordinates++) = FT(0);
        return std::make_pair(coordinates, true);
      }

      case Query_point_location::ON_EDGE: {
        CGAL_assertion(index >= 0 && index < n);

        if (index == n - 1)
          return coordinates_on_last_edge_2(
            polygon, query, coordinates, traits, vertex_map);

        const std::size_t indexp = (index + 1) % n;
        const auto& source = get(vertex_map, *(polygon.begin() + index));
        const auto& target = get(vertex_map, *(polygon.begin() + indexp));

        for (std::size_t i = 0; i < n; ++i)
          if (i == index) {
            linear_coordinates_2(source, target, query, coordinates, traits); ++i;
          } else {
            *(coordinates++) = FT(0);
          }
        return std::make_pair(coordinates, true);
      }

      default: {
        internal::get_default(n, coordinates);
        return std::make_pair(coordinates, false);
      }
    }
    return std::make_pair(coordinates, false);
  }

  // Do we need it at all? Can we use DH weights inside Harmonic coordinates?
  template<typename GeomTraits>
  typename GeomTraits::FT cotangent_2(
    const typename GeomTraits::Vector_2& v1,
    const typename GeomTraits::Vector_2& v2,
    const GeomTraits traits) {

    using FT = typename GeomTraits::FT;
    const auto cross_product_2 = traits.compute_determinant_2_object();
    const auto scalar_product_2 = traits.compute_scalar_product_2_object();

    const FT det = cross_product_2(v1, v2);
    const FT dot = scalar_product_2(v1, v2);

    const FT cot_denominator = CGAL::abs(det);
    CGAL_assertion(cot_denominator != FT(0));
    return dot / cot_denominator;
  }

  // This function is taken from the Polygon_2_algorithms.h header.
  // But it is updated to support property maps.
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap>
  bool is_simple_2(
    const Polygon& polygon,
    const GeomTraits traits,
    const VertexMap vertex_map) {

    const auto first = polygon.begin();
    const auto last = polygon.end();
    if (first == last) return true;

    std::vector<typename GeomTraits::Point_2> poly;
    poly.reserve(polygon.size());
    for (const auto& vertex : polygon)
      poly.push_back(get(vertex_map, vertex));
    return CGAL::is_simple_2(poly.begin(), poly.end(), traits);
  }

} // namespace internal
} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_UTILS_2_H
