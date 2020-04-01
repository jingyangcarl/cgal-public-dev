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

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/assertions.h>

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename FT>
  static FT tolerance() {
    return FT(1) / FT(100000);
  }

  template<typename FT>
  static FT max_value() {
    return FT(1000000000000);
  }

  template<typename Point_2>
  Point_2 middle_point_2(
    const Point_2& source, const Point_2& target) {
    
    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    const FT half = FT(1) / FT(2);
    const FT x = half * (source.x() + target.x());
    const FT y = half * (source.y() + target.y());
    return Point_2(x, y);
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::FT
  length_2(const Segment_2& segment) { 
    
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    
    return static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(segment.squared_length())));
  }

  template<typename Vector>
  void normalize(Vector& v) {
    
    using Traits = typename Kernel_traits<Vector>::Kernel;
    using FT = typename Traits::FT;
    
    v /= static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(v.squared_length())));
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::Vector_2
  direction_2(const Segment_2& segment) { 
    
    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;
    using Vector_2 = typename Traits::Vector_2;

    Vector_2 v = segment.to_vector();
    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0))) 
      v = -v;
    normalize(v);
    return v;
  }

  template<typename Vector_2>
  typename Kernel_traits<Vector_2>::Kernel::FT
  orientation_2(const Vector_2& v) {
    
    using Traits = typename Kernel_traits<Vector_2>::Kernel;
    using FT = typename Traits::FT;

    const FT atan_rad = static_cast<FT>(
      std::atan2(CGAL::to_double(v.y()), CGAL::to_double(v.x())));
    FT angle_deg = atan_rad * FT(180) / static_cast<FT>(CGAL_PI);
    if (angle_deg < FT(0)) 
      angle_deg += FT(180);
    return angle_deg;
  }

  template<
  typename Point_2, 
  typename FT>
  Point_2 transform_coordinates_2(
    const Point_2& barycenter, 
    const Point_2& frame_origin, 
    const FT angle_deg) {

    const double angle_rad = 
      CGAL_PI * CGAL::to_double(angle_deg) / 180.0;

    const FT cos_val = static_cast<FT>(std::cos(angle_rad));
    const FT sin_val = static_cast<FT>(std::sin(angle_rad));

    const FT diff_x = barycenter.x() - frame_origin.x();
    const FT diff_y = barycenter.y() - frame_origin.y();

    const FT x = diff_x * cos_val + diff_y * sin_val;
    const FT y = diff_y * cos_val - diff_x * sin_val;

    return Point_2(x, y);
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::FT
  compute_angle_2(
    const Segment_2& longest, 
    const Segment_2& segment) {

    using Traits = typename Kernel_traits<Segment_2>::Kernel;
    using FT = typename Traits::FT;

    const auto v1 =  segment.to_vector();
    const auto v2 = -longest.to_vector();

    const FT det = CGAL::determinant(v1, v2);
    const FT dot = CGAL::scalar_product(v1, v2);
    const FT angle_rad = static_cast<FT>(
      std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
    const FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI); 
    return angle_deg;
  }

  template<typename FT>
  FT convert_angle_2(const FT angle) {
    
    FT angle_2 = angle;
    if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
    else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
    return angle_2;
  }

  template<typename Segment_2>
  typename Kernel_traits<Segment_2>::Kernel::FT
  angle_2_degrees(
    const Segment_2& longest,
    const Segment_2& segment) {

    const auto angle_2 = compute_angle_2(
      longest, segment);
    return convert_angle_2(angle_2);
  }

  template<
  typename FT,
  typename Point_2>
  void rotate_point_2(
    const FT angle, 
    const Point_2& barycenter, 
    Point_2& p) {

		FT x = p.x(); x -= barycenter.x();
		FT y = p.y(); y -= barycenter.y();
		
    p = Point_2(x, y);
    const double tmp_angle = CGAL::to_double(angle);
    const FT c = static_cast<FT>(std::cos(tmp_angle));
		const FT s = static_cast<FT>(std::sin(tmp_angle));

		x = p.x() * c - p.y() * s; x += barycenter.x();
		y = p.y() * c + p.x() * s; y += barycenter.y();
		p = Point_2(x, y);
	} 

  template<typename Point_2>
  Point_2 barycenter_2(
    const std::vector<Point_2>& points) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    using FT = typename Traits::FT;

    CGAL_assertion(points.size() > 0);
    FT x = FT(0), y = FT(0);
    for (const auto& p : points) {
      x += p.x();
      y += p.y();
    }
    x /= static_cast<FT>(points.size());
    y /= static_cast<FT>(points.size());
    return Point_2(x, y);
  }

} // internal
} // Shape_regularization
} // CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
