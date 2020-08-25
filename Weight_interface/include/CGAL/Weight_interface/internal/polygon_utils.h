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

#ifndef CGAL_WEIGHT_INTERFACE_INTERNAL_POLYGON_UTILS_H
#define CGAL_WEIGHT_INTERFACE_INTERNAL_POLYGON_UTILS_H

// #include <CGAL/license/Weight_interface.h>

// STL includes.
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace CGAL {
namespace Weights {
namespace internal {

  // Normalize values.
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

  // Polygon type enum.
  enum class Polygon_type {

    // Concave polygon = non-convex polygon.
    CONCAVE = 0,

    // This is a convex polygon with collinear vertices.
    WEAKLY_CONVEX = 1,

    // This is a convex polygon without collinear vertices.
    STRICTLY_CONVEX = 2
  };

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

      // Test all the consequent triplets of polygon vertices on collinearity.
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
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_INTERNAL_POLYGON_UTILS_H
