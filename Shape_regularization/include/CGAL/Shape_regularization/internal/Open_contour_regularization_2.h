// Copyright (c) 2019 GeometryFactory Sarl (France).
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

#ifndef CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Shape_regularization/enum.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Open_contour_regularization_2 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    using FT_pair = std::pair<FT, FT>;

    Open_contour_regularization_2(
      Input_range& input_range,
      Point_map point_map) :
    m_input_range(input_range),
    m_point_map(point_map) { 
      
      CGAL_precondition(input_range.size() >= 2);
    }

    void set_principal_directions(
      const std::vector<std::size_t>& directions) {

      CGAL_precondition(
        directions.size() == m_input_range.size() - 1);
      m_group = directions;
    }

    void estimate_principal_directions(
      const Direction_type direction_type,
      const FT min_length_2,
      const FT max_angle_2) {

      switch (direction_type) {
        case Direction_type::LONGEST:
          set_longest_direction();
          break;
        case Direction_type::MULTIPLE:
          set_multiple_directions(
            min_length_2, max_angle_2);
          break;
        default:
          set_longest_direction();
          break;
      };
    }

    void regularize() {
      
    }

    std::size_t number_of_principal_directions() const {

      std::set<std::size_t> unique;
      for (std::size_t index : m_group)
        if (index != std::size_t(-1))
          unique.insert(index);
      return unique.size();
    }

  private:
    Input_range& m_input_range;
    Point_map m_point_map;

    std::vector<FT_pair> m_bounds;
    std::vector<Segment_2> m_longest;
    std::vector<std::size_t> m_group;

    void set_longest_direction() {

    }

    void set_multiple_directions(
      const FT min_length_2,
      const FT max_angle_2) {
      
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OPEN_CONTOUR_REGULARIZATION_2_H
