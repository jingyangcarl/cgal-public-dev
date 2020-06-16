// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_GENERALIZED_WEIGHTS_UTILS_2_H
#define CGAL_GENERALIZED_WEIGHTS_UTILS_2_H

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

namespace CGAL {
namespace Generalized_weights {
namespace internal {

// Raises value to the power.
template<typename GeomTraits>
const typename GeomTraits::FT power(
  const GeomTraits& traits,
  const typename GeomTraits::FT value,
  const typename GeomTraits::FT p) {

  using FT = typename GeomTraits::FT;
  const double base = CGAL::to_double(value);
  const double exp = CGAL::to_double(p);
  return static_cast<FT>(std::pow(base, exp));
}

// Computes distance between two 2D points.
template<typename GeomTraits>
const typename GeomTraits::FT distance_2(
  const GeomTraits& traits,
  const typename GeomTraits::Point_2& p,
  const typename GeomTraits::Point_2& q) {

  using FT = typename GeomTraits::FT;
  const auto squared_distance_2 =
    traits.compute_squared_distance_2_object();
  return static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_distance_2(p, q))));
}

// Computes length of a 2D vector.
template<typename GeomTraits>
const typename GeomTraits::FT length_2(
  const GeomTraits& traits,
  const typename GeomTraits::Vector_2& v) {

  using FT = typename GeomTraits::FT;
  const auto squared_length_2 =
    traits.compute_squared_length_2_object();
  return static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_length_2(v))));
}

// Normalizes a 2D vector.
template<typename GeomTraits>
void normalize_2(
  const GeomTraits& traits,
  typename GeomTraits::Vector_2& v) {

  using FT = typename GeomTraits::FT;
  const FT length = length_2(traits, v);
  CGAL_assertion(length != FT(0));
  if (length == FT(0)) return;
  v /= length;
}

// Computes cotanget between two 2D vectors.
template<typename GeomTraits>
const typename GeomTraits::FT cotangent_2(
  const GeomTraits& traits,
  const typename GeomTraits::Point_2& p,
  const typename GeomTraits::Point_2& q,
  const typename GeomTraits::Point_2& r) {

  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  const auto dot_product_2 =
    traits.compute_scalar_product_2_object();
  const auto cross_product_2 =
    traits.compute_determinant_2_object();

  const Vector_2 u = Vector_2(q, r);
  const Vector_2 v = Vector_2(q, p);

  const FT dot = dot_product_2(u, v);
  const FT cross = cross_product_2(u, v);

  const FT length = CGAL::abs(cross);
  CGAL_assertion(length != FT(0));
  if (length != FT(0)) return dot / length;
  else return FT(0); // undefined
}

// Computes tanget between two 2D vectors.
template<typename GeomTraits>
const typename GeomTraits::FT tangent_2(
  const GeomTraits& traits,
  const typename GeomTraits::Point_2& p,
  const typename GeomTraits::Point_2& q,
  const typename GeomTraits::Point_2& r) {

  using FT = typename GeomTraits::FT;
  using Vector_2 = typename GeomTraits::Vector_2;

  const auto dot_product_2 =
    traits.compute_scalar_product_2_object();
  const auto cross_product_2 =
    traits.compute_determinant_2_object();

  const Vector_2 u = Vector_2(q, r);
  const Vector_2 v = Vector_2(q, p);

  const FT dot = dot_product_2(u, v);
  const FT cross = cross_product_2(u, v);

  const FT length = CGAL::abs(cross);
  CGAL_assertion(dot != FT(0));
  if (dot != FT(0)) return length / dot;
  else return FT(0); // undefined
}

// Computes distance between two 3D points.
template<typename GeomTraits>
const typename GeomTraits::FT distance_3(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& p,
  const typename GeomTraits::Point_3& q) {

  using FT = typename GeomTraits::FT;
  const auto squared_distance_3 =
    traits.compute_squared_distance_3_object();
  return static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_distance_3(p, q))));
}

// Computes length of a 3D vector.
template<typename GeomTraits>
const typename GeomTraits::FT length_3(
  const GeomTraits& traits,
  const typename GeomTraits::Vector_3& v) {

  using FT = typename GeomTraits::FT;
  const auto squared_length_3 =
    traits.compute_squared_length_3_object();
  return static_cast<FT>(
    CGAL::sqrt(CGAL::to_double(squared_length_3(v))));
}

// Normalizes a 3D vector.
template<typename GeomTraits>
void normalize_3(
  const GeomTraits& traits,
  typename GeomTraits::Vector_3& v) {

  using FT = typename GeomTraits::FT;
  const FT length = length_3(traits, v);
  CGAL_assertion(length != FT(0));
  if (length == FT(0)) return;
  v /= length;
}

// Computes cotanget between two 3D vectors.
template<typename GeomTraits>
const typename GeomTraits::FT cotangent_3(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& p,
  const typename GeomTraits::Point_3& q,
  const typename GeomTraits::Point_3& r) {

  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();
  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();

  const Vector_3 u = Vector_3(q, r);
  const Vector_3 v = Vector_3(q, p);

  const FT dot = dot_product_3(u, v);
  const Vector_3 cross = cross_product_3(u, v);

  const FT length = length_3(traits, cross);
  CGAL_assertion(length != FT(0));
  if (length != FT(0)) return dot / length;
  else return FT(0); // undefined
}

// Computes tanget between two 3D vectors.
template<typename GeomTraits>
const typename GeomTraits::FT tangent_3(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& p,
  const typename GeomTraits::Point_3& q,
  const typename GeomTraits::Point_3& r) {

  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();
  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();

  const Vector_3 u = Vector_3(q, r);
  const Vector_3 v = Vector_3(q, p);

  const FT dot = dot_product_3(u, v);
  const Vector_3 cross = cross_product_3(u, v);

  const FT length = length_3(traits, cross);
  CGAL_assertion(dot != FT(0));
  if (dot != FT(0)) return length / dot;
  else return FT(0); // undefined
}

// Computes area of a 3D triangle.
template<typename GeomTraits>
const typename GeomTraits::FT area_3(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& p,
  const typename GeomTraits::Point_3& q,
  const typename GeomTraits::Point_3& r) {

  using FT = typename GeomTraits::FT;
  using Vector_3 = typename GeomTraits::Vector_3;

  const Vector_3 u = Vector_3(q, r);
  const Vector_3 v = Vector_3(q, p);

  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();
  const Vector_3 cross = cross_product_3(u, v);
  const FT half = FT(1) / FT(2);
  const FT area = half * length_3(traits, cross);
  return area;
}

// Computes 3D barycenter of the points a, b, c, and d.
template<typename GeomTraits>
void barycenter_3(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& a,
  const typename GeomTraits::Point_3& b,
  const typename GeomTraits::Point_3& c,
  const typename GeomTraits::Point_3& d,
  typename GeomTraits::Point_3& barycenter) {

  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  const FT quater = FT(1) / FT(4);
  const FT x = quater * (a.x() + b.x() + c.x() + d.x());
  const FT y = quater * (a.y() + b.y() + c.y() + d.y());
  const FT z = quater * (a.z() + b.z() + c.z() + d.z());
  barycenter = Point_3(x, y, z);
}

// Computes 3D angle between two vectors.
template<typename GeomTraits>
const double angle_3(
  const GeomTraits& traits,
  const typename GeomTraits::Vector_3& v1,
  const typename GeomTraits::Vector_3& v2) {

  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();
  const double dot =
    CGAL::to_double(dot_product_3(v1, v2));

  double angle_rad = 0.0;
  if (dot < -1.0) angle_rad = std::acos(-1.0);
  else if (dot > 1.0) angle_rad = std::acos(1.0);
  angle_rad = std::acos(dot);
  return angle_rad;
}

// Rotates a 3D point around axis.
template<typename GeomTraits>
const typename GeomTraits::Point_3 rotate_point_3(
  const GeomTraits& traits,
  const double angle_rad,
  const typename GeomTraits::Vector_3& axis,
  const typename GeomTraits::Point_3& query) {

  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;

  const FT c = static_cast<FT>(std::cos(angle_rad));
  const FT s = static_cast<FT>(std::sin(angle_rad));
  const FT C = FT(1) - c;

  const auto x = axis.x();
  const auto y = axis.y();
  const auto z = axis.z();

  return Point_3(
    (x * x * C + c)     * query.x() +
    (x * y * C - z * s) * query.y() +
    (x * z * C + y * s) * query.z(),
    (y * x * C + z * s) * query.x() +
    (y * y * C + c)     * query.y() +
    (y * z * C - x * s) * query.z(),
    (z * x * C - y * s) * query.x() +
    (z * y * C + x * s) * query.y() +
    (z * z * C + c)     * query.z());
}

// Computes two 3D orthogonal base vectors wrt a given normal.
template<typename GeomTraits>
void orthogonal_bases_3(
  const GeomTraits& traits,
  const typename GeomTraits::Vector_3& normal,
  typename GeomTraits::Vector_3& b1,
  typename GeomTraits::Vector_3& b2) {

  using Vector_3 = typename GeomTraits::Vector_3;
  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();

  const auto nx = normal.x();
  const auto ny = normal.y();
  const auto nz = normal.z();

  if (CGAL::abs(nz) >= CGAL::abs(ny))
    b1 = Vector_3(nz, 0, -nx);
  else
    b1 = Vector_3(ny, -nx, 0);
  b2 = cross_product_3(normal, b1);

  normalize_3(traits, b1);
  normalize_3(traits, b2);
}

// Converts a 3D point into a 2D point wrt to a given plane.
template<typename GeomTraits>
const typename GeomTraits::Point_2 to_2d(
  const GeomTraits& traits,
  const typename GeomTraits::Vector_3& b1,
  const typename GeomTraits::Vector_3& b2,
  const typename GeomTraits::Point_3& origin,
  const typename GeomTraits::Point_3& query) {

  using Point_2  = typename GeomTraits::Point_2;
  using Vector_3 = typename GeomTraits::Vector_3;
  const auto dot_product_3 =
    traits.compute_scalar_product_3_object();

  const Vector_3 v = Vector_3(origin, query);
  const auto x = dot_product_3(b1, v);
  const auto y = dot_product_3(b2, v);
  return Point_2(x, y);
}

// Flattens an arbitrary quad into a planar quad.
template<typename GeomTraits>
void flatten(
  const GeomTraits& traits,
  const typename GeomTraits::Point_3& query, // query point
  const typename GeomTraits::Point_3& nm,    // prev neighbor/vertex
  const typename GeomTraits::Point_3& nj,    // curr neighbor/vertex
  const typename GeomTraits::Point_3& np,    // next neighbor/vertex
  typename GeomTraits::Point_2& pq,
  typename GeomTraits::Point_2& pm,
  typename GeomTraits::Point_2& pj,
  typename GeomTraits::Point_2& pp) {

  // std::cout << std::endl;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  const auto cross_product_3 =
    traits.construct_cross_product_vector_3_object();

  // Compute barycenter.
  Point_3 center;
  barycenter_3(traits, query, nm, nj, np, center);
  // std::cout << "barycenter: " << center << std::endl;

  // Translate.
  const Point_3 vq = Point_3(
    query.x() - center.x(), query.y() - center.y(), query.z() - center.z());
  const Point_3 vm = Point_3(
    nm.x() - center.x(), nm.y() - center.y(), nm.z() - center.z());
  const Point_3 vj = Point_3(
    nj.x() - center.x(), nj.y() - center.y(), nj.z() - center.z());
  const Point_3 vp = Point_3(
    np.x() - center.x(), np.y() - center.y(), np.z() - center.z());

  // std::cout << "tr vq: " << vq << std::endl;
  // std::cout << "tr vm: " << vm << std::endl;
  // std::cout << "tr vj: " << vj << std::endl;
  // std::cout << "tr vp: " << vp << std::endl;

  // Middle axis.
  Vector_3 ax = Vector_3(vq, vj);
  normalize_3(traits, ax);

  // Prev and next vectors.
  Vector_3 v1 = Vector_3(vq, vm);
  Vector_3 v2 = Vector_3(vq, vp);

  normalize_3(traits, v1);
  normalize_3(traits, v2);

  // Two triangle normals.
  Vector_3 n1 = cross_product_3(v1, ax);
  Vector_3 n2 = cross_product_3(ax, v2);

  normalize_3(traits, n1);
  normalize_3(traits, n2);

  // std::cout << "normal n1: " << n1 << std::endl;
  // std::cout << "normal n2: " << n2 << std::endl;

  // Angle between two normals.
  const double angle_rad = angle_3(traits, n1, n2);
  // std::cout << "angle deg n1 <-> n2: " << angle_rad * 180.0 / CGAL_PI << std::endl;

  // Rotate vp around ax so that it lands onto the plane [query, vm, vj].
  const Point_3& rq = vq;
  const Point_3& rm = vm;
  const Point_3& rj = vj;
  const Point_3  rp = rotate_point_3(traits, angle_rad, ax, vp);
  // std::cout << "rotated vp: " << rp << std::endl;

  // Compute orthogonal base vectors.
  Vector_3 b1, b2;
  const auto& normal = n1;
  orthogonal_bases_3(traits, normal, b1, b2);

  // const auto angle12 = angle_3(traits, b1, b2);
  // std::cout << "angle deg b1 <-> b2: " << angle12 * 180.0 / CGAL_PI << std::endl;

  // Flatten a quad.
  const auto& origin = vq;
  pq = to_2d(traits, b1, b2, origin, rq);
  pm = to_2d(traits, b1, b2, origin, rm);
  pj = to_2d(traits, b1, b2, origin, rj);
  pp = to_2d(traits, b1, b2, origin, rp);

  // std::cout << "vq: " << pq << std::endl;
  // std::cout << "vm: " << pm << std::endl;
  // std::cout << "vj: " << pj << std::endl;
  // std::cout << "vp: " << pp << std::endl;

  // const auto area_2 = traits.compute_area_2_object();
  // std::cout << "Am: " << area_2(pq, pm, pj) << std::endl;
  // std::cout << "Aj: " << area_2(pq, pj, pp) << std::endl;
  // std::cout << "C: "  << area_2(pm, pj, pp) << std::endl;
  // std::cout << "B: "  << area_2(pq, pm, pp) << std::endl;
}

} // namespace internal
} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_UTILS_2_H
