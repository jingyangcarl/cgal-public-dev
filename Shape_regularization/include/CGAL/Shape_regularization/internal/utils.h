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
// Author(s)     : Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <tuple>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Aff_transformation_2.h>
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

  template<
  typename FT,
  typename Point_2>
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

  template<typename Vector_d>
  void normalize(Vector_d& v) {

    using Traits = typename Kernel_traits<Vector_d>::Kernel;
    using FT = typename Traits::FT;

    v /= static_cast<FT>(
      CGAL::sqrt(CGAL::to_double(v.squared_length())));
  }

  template<typename Vector_2>
  typename Kernel_traits<Vector_2>::Kernel::Direction_2
  direction_2(Vector_2& v) {

    using Traits = typename Kernel_traits<Vector_2>::Kernel;
    using FT = typename Traits::FT;
    using Direction_2 = typename Traits::Direction_2;

    if (v.y() < FT(0) || (v.y() == FT(0) && v.x() < FT(0)))
      v = -v;
    normalize(v);
    return Direction_2(v);
  }

  template<typename Vector_2>
  typename Kernel_traits<Vector_2>::Kernel::FT
  orientation_2(const Vector_2& v) {

    using Traits = typename Kernel_traits<Vector_2>::Kernel;
    using FT = typename Traits::FT;

    const FT angle_rad = static_cast<FT>(std::atan2(
      CGAL::to_double(v.y()),
      CGAL::to_double(v.x())));

    FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
    if (angle_deg < FT(0))
      angle_deg += FT(180);
    return angle_deg;
  }

  template<typename FT>
  FT orientation_to_angle_2(
    const FT oi, const FT oj) {

    const FT diff_ij = oi - oj;
    const double diff_90 = std::floor(CGAL::to_double(diff_ij / FT(90)));
    const FT to_lower = FT(90) *  static_cast<FT>(diff_90)          - diff_ij;
    const FT to_upper = FT(90) * (static_cast<FT>(diff_90) + FT(1)) - diff_ij;

    const FT angle_deg =
      CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;
    return angle_deg;
  }

  // It is used only in the contour regularization.
  template<typename Direction_2>
  typename Kernel_traits<Direction_2>::Kernel::FT
  invar90_angle_2(
    const Direction_2& di,
    const Direction_2& dj) {

    using Traits = typename Kernel_traits<Direction_2>::Kernel;
    using FT = typename Traits::FT;

    const auto vdi = di.to_vector();
    const FT oi = orientation_2(vdi);
    const auto vdj = dj.to_vector();
    const FT oj = orientation_2(vdj);

    return orientation_to_angle_2(oi, oj);
  }

  template<typename Direction_2>
  typename Kernel_traits<Direction_2>::Kernel::FT
  compute_angle_2(
    const Direction_2& reference,
    const Direction_2& direction) {

    using Traits = typename Kernel_traits<Direction_2>::Kernel;
    using FT = typename Traits::FT;

    const auto v1 =  direction.to_vector();
    const auto v2 = -reference.to_vector();

    const FT det = CGAL::determinant(v1, v2);
    const FT dot = CGAL::scalar_product(v1, v2);
    const FT angle_rad = static_cast<FT>(
      std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
    const FT angle_deg = angle_rad * FT(180) / static_cast<FT>(CGAL_PI);
    return angle_deg;
  }

  template<typename FT>
  FT convert_angle_2(const FT angle_2) {

    FT angle = angle_2;
    if (angle > FT(90)) angle = FT(180) - angle;
    else if (angle < -FT(90)) angle = FT(180) + angle;
    return angle;
  }

  template<typename Direction_2>
  typename Kernel_traits<Direction_2>::Kernel::FT
  angle_2(
    const Direction_2& reference,
    const Direction_2& direction) {

    const auto angle_2 = compute_angle_2(
      reference, direction);
    return convert_angle_2(angle_2);
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

  // Redo this function via CGAL affine transform!
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

  // Redo this function via CGAL affine transform!
  template<
  typename FT,
  typename Segment_2>
  void rotate_segment_2(
    const FT angle_2_deg,
    const FT ref_angle_2_deg,
    Segment_2& segment) {

    FT angle_deg = angle_2_deg;
    if (angle_deg < FT(0)) angle_deg += ref_angle_2_deg;
    else if (angle_deg > FT(0)) angle_deg -= ref_angle_2_deg;

    auto source = segment.source();
    auto target = segment.target();

    const auto barycenter = internal::middle_point_2(source, target);
    const FT angle_rad = angle_deg * static_cast<FT>(CGAL_PI) / FT(180);

    rotate_point_2(angle_rad, barycenter, source);
    rotate_point_2(angle_rad, barycenter, target);
    segment = Segment_2(source, target);
  }

  template<
  typename FT,
  typename Direction_2>
  void rotate_direction_2(
    const FT angle_deg,
    Direction_2& direction) {

    using Traits = typename Kernel_traits<Direction_2>::Kernel;
    using Transformation_2 = typename Traits::Aff_transformation_2;

    const FT angle_rad = angle_deg * static_cast<FT>(CGAL_PI) / FT(180);
    const double sinval = std::sin(CGAL::to_double(angle_rad));
    const double cosval = std::cos(CGAL::to_double(angle_rad));
    const Transformation_2 rotate_2(CGAL::ROTATION, sinval, cosval);
    direction = rotate_2(direction);
  }

  template<typename FT>
  double radians_2(FT angle) {

    if (angle < FT(0)) angle += FT(180);
    else if (angle > FT(180)) angle -= FT(180);
    const double angle_rad = CGAL::to_double(
      angle * static_cast<FT>(CGAL_PI) / FT(180));
    return angle_rad;
  }

  // Remove this function since max_angle cannot be equal 0!
  template<typename Segment_2>
  std::size_t key_angle_2(
    const double max_angle,
    const Segment_2& segment) {

    auto v = segment.to_vector();
    const auto direction = internal::direction_2(v).to_vector();
    const auto orientation = internal::orientation_2(direction);
    double fvalue = std::ceil(CGAL::to_double(orientation));
    if (fvalue >= 180.0) fvalue -= 180.0;
    CGAL_assertion(max_angle != 0.0);
    const std::size_t num = static_cast<std::size_t>(
      std::floor(fvalue / max_angle));
    const std::size_t angle = static_cast<std::size_t>(
      num * max_angle);
    return angle;
  }

} // internal
} // Shape_regularization
} // CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
