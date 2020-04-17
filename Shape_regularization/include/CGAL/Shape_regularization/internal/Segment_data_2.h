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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2_H

// #include <CGAL/license/Shape_regularization.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  struct Segment_data_2 {

  public:
    using Traits = GeomTraits;
    using Segment_2 = typename Traits::Segment_2;
    using Point_2 = typename Traits::Point_2;
    using Vector_2  = typename Traits::Vector_2;
    using FT = typename Traits::FT;

    const Segment_2& segment;
    const std::size_t index;
    
    Point_2  barycenter;
    FT       length;
    Vector_2 direction;
    FT       orientation;
    Point_2  ref_coords;
    FT       a, b, c;

    Segment_data_2(
      const Segment_2& s,
      const std::size_t sindex):
    segment(s), index(sindex) {

      barycenter = internal::middle_point_2(
        segment.source(), segment.target());
      length = internal::length_2(segment);
      direction = internal::direction_2(segment);
      orientation = internal::orientation_2(direction);

      const double angle_rad = CGAL::to_double(
        orientation * static_cast<FT>(CGAL_PI) / FT(180));
      a = -static_cast<FT>(std::sin(angle_rad));
      b =  static_cast<FT>(std::cos(angle_rad));
      c = -a * barycenter.x() - b * barycenter.y();
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_SEGMENT_DATA_2_H
