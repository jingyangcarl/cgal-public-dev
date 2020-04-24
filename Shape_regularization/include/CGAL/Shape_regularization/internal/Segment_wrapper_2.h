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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  struct Segment_wrapper_2 {

  public:
    using Traits = GeomTraits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2  = typename Traits::Direction_2;
    
    std::size_t index = std::size_t(-1);
    Point_2 barycenter, ref_coords;
    FT orientation, length, a, b, c;
    Direction_2 direction;
    bool is_used = false;

    void set_index(const std::size_t index_) {
      index = index_; is_used = true;
    }

    void set_length(
      const Segment_2& segment) {
      length = internal::length_2(segment);
    }

    void set_barycenter(
      const Segment_2& segment) {

      barycenter = internal::middle_point_2(
        segment.source(), segment.target());
    }

    void set_direction(
      const Segment_2& segment) {
      direction = internal::direction_2(segment); 
    }

    void set_orientation(
      const Segment_2& segment) {

      set_direction(segment);
      const auto dir = direction.to_vector();
      orientation = internal::orientation_2(dir);
    }

    void set_abc(
      const Segment_2& segment) {
      
      set_barycenter(segment);
      set_orientation(segment);
      
      const double angle_rad = CGAL::to_double(
        orientation * static_cast<FT>(CGAL_PI) / FT(180));
      a = -static_cast<FT>(std::sin(angle_rad));
      b = +static_cast<FT>(std::cos(angle_rad));
      c = -a * barycenter.x() - b * barycenter.y();
    }

    void set_all(
      const std::size_t index_,
      const Segment_2& segment) {
      
      set_index(index_);
      set_length(segment);
      set_abc(segment);
    }

    void set_ref_coords(
      const Point_2& frame_origin) {

      ref_coords = internal::transform_coordinates_2(
        barycenter, frame_origin, orientation);
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_WRAPPER_2_H
