// Copyright (c) 2020 GeometryFactory Sarl (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Andreas Fabri
//

#ifndef CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H
#define CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Longest_principal_direction_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;

    using FT_pair = std::pair<FT, FT>;

    Longest_principal_direction_2(
      const InputRange& input_range,
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_point_map(point_map)
    { }

    void estimate(
      const bool is_closed,
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      if (is_closed)
        estimate_closed(bounds, directions, assigned);
      else 
        estimate_open(bounds, directions, assigned);
    }

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;

    void estimate_closed(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction_closed();

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(m_input_range.size(), 0);
    }

    Direction_2 compute_longest_direction_closed() const {

      FT max_length = -FT(1);
      std::size_t index = std::size_t(-1);

      CGAL_assertion(m_input_range.size() >= 3);
      const std::size_t n = m_input_range.size();
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;

        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        const FT sq_length = CGAL::squared_distance(source, target);
        if (sq_length > max_length) {
          index = i; max_length = sq_length;
        }
      }

      const std::size_t i = index;
      const std::size_t ip = (i + 1) % n;
      const auto& source = get(m_point_map, *(m_input_range.begin() + i));
      const auto& target = get(m_point_map, *(m_input_range.begin() + ip));
      
      const Segment_2 segment = Segment_2(source, target);
      return internal::segment_to_direction_2(segment);
    }

    void estimate_open(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      CGAL_assertion(m_input_range.size() >= 2);

      bounds.clear(); bounds.resize(1);
      bounds[0] = std::make_pair(FT(45), FT(45));

      directions.clear(); directions.resize(1);
      directions[0] = compute_longest_direction_open();

      // 0 is the index of the direction in the `directions`.
      assigned.clear();
      assigned.resize(m_input_range.size() - 1, 0);
    }

    Direction_2 compute_longest_direction_open() const {

      FT max_length = -FT(1);
      std::size_t index = std::size_t(-1);

      CGAL_assertion(m_input_range.size() >= 2);
      const std::size_t n = m_input_range.size();
      for (std::size_t i = 0; i < n - 1; ++i) {
        const std::size_t ip = i + 1;

        const auto& source = get(m_point_map, *(m_input_range.begin() + i));
        const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

        const FT sq_length = CGAL::squared_distance(source, target);
        if (sq_length > max_length) {
          index = i; max_length = sq_length;
        }
      }

      const std::size_t i = index;
      const std::size_t ip = i + 1;
      const auto& source = get(m_point_map, *(m_input_range.begin() + i));
      const auto& target = get(m_point_map, *(m_input_range.begin() + ip));

      const Segment_2 segment = Segment_2(source, target);
      return internal::segment_to_direction_2(segment);
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_LONGEST_PRINCIPAL_DIRECTION_2_H
