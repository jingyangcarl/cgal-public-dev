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

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_2.h>

// Internal includes.
#include <CGAL/Shape_regularization/enum.h>
#include <CGAL/Shape_regularization/internal/utils.h>

// Saver.
#include "../../../../examples/Shape_regularization/include/Saver.h"

// TODO:
// * Can I further simplify this class?
// * Can I use squared distance here instead of distance?
// * Improve find_central_segment().
// * Improve orth segments, which are added during an optimization step. They 
// are too far away from the correct place since they are placed in the middle of
// the segment.
// * Do we actually need make_segments_collienar() in the contour connection?
// * Improve intersection.

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class Contour_regularization_base_2 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    struct Segment_wrapper_2 {

      Segment_2 segment;
      std::size_t index = std::size_t(-1);
      std::size_t group = std::size_t(-1);
      bool is_valid_direction = false;
      bool is_used = false;
    };

    using FT_pair = std::pair<FT, FT>;
    using Segments_2 = std::vector<Segment_2>;
    using Segment_wrappers_2 = std::vector<Segment_wrapper_2>;

    Contour_regularization_base_2(
      const FT min_length_2,
      const FT max_angle_2,
      const FT max_offset_2) :
    m_min_length_2(min_length_2),
    m_max_angle_2(max_angle_2),
    m_max_offset_2(max_offset_2),
    m_angle_threshold_2(FT(5)),
    m_verbose(true) { 

      CGAL_precondition(min_length_2 >= FT(0));
      CGAL_precondition(max_angle_2 >= FT(0) && max_angle_2 <= FT(90));
      CGAL_precondition(max_offset_2 >= FT(0));
    }

    const FT get_min_length_2() const {
      return m_min_length_2;
    }

    const FT get_max_angle_2() const {
      return m_max_angle_2;
    }

    const FT get_max_offset_2() const {
      return m_max_offset_2;
    }

    const FT get_angle_threshold_2() const {
      return m_angle_threshold_2;
    }

    const bool verbose() const {
      return m_verbose;
    }

    const bool is_valid_principal_direction(
      const Segment_2& segment) const {
      
      const FT threshold = get_min_length_2() * FT(2);
      const FT squared_threshold = threshold * threshold;
      return segment.squared_length() >= squared_threshold;
    }

    Direction_2 compute_longest_direction(
      const std::vector<Segment_wrapper_2>& wraps) const {

      const std::size_t n = wraps.size();
      CGAL_assertion(n != 0);

      FT max_length = -FT(1);
      std::size_t longest = std::size_t(-1);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& wrap = wraps[i];
        const FT sq_length = wrap.segment.squared_length();
        if (sq_length > max_length) {
          longest = i; max_length = sq_length;
        }
      }

      CGAL_assertion(longest != std::size_t(-1));
      return Direction_2(wraps[longest].segment);
    }

  private:
    const FT m_min_length_2;
    const FT m_max_angle_2;
    const FT m_max_offset_2;

    const FT m_angle_threshold_2;
    const bool m_verbose;
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_BASE_2_H
