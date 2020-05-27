// Copyright (c) 2015-2020 INRIA Sophia-Antipolis and GeometryFactory SARL (France).
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
// Author(s)     : Florent Lafarge, Simon Giraudot, Gennadii Sytov, Dmitry Anisimov
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
#include <CGAL/centroid.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/squared_distance_3.h>
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

  template<typename Traits>
  struct Plane_cluster {

    bool is_free;
    std::vector<std::size_t> planes;
    std::vector<std::size_t> coplanar_group;
    std::vector<std::size_t> orthogonal_clusters;
    typename Traits::Vector_3 normal;
    typename Traits::FT cosangle_symmetry;
    typename Traits::FT area;
    typename Traits::FT cosangle_centroid;

    Plane_cluster():
      is_free(true),
      normal(
        typename Traits::FT(0),
        typename Traits::FT(0),
        typename Traits::FT(1)),
      cosangle_symmetry(typename Traits::FT(0)),
      area(typename Traits::FT(0)),
      cosangle_centroid(typename Traits::FT(0))
    { }
  };

  template<typename Traits>
  typename Traits::Vector_3 regularize_normal(
    const typename Traits::Vector_3& n,
    const typename Traits::Vector_3& symmetry_direction,
    const typename Traits::FT cos_symmetry) {

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Line_3 Line;
    typedef typename Traits::Plane_3 Plane;

    const Point pt_symmetry = CGAL::ORIGIN + cos_symmetry * symmetry_direction;
    const Plane plane_symmetry(pt_symmetry, symmetry_direction);
    const Point pt_normal = CGAL::ORIGIN + n;

    if (n != symmetry_direction || n != -symmetry_direction) {
      const Plane plane_cut(
        CGAL::ORIGIN, pt_normal, CGAL::ORIGIN + symmetry_direction);

      Line line;
      const CGAL::Object ob_1 = CGAL::intersection(plane_cut, plane_symmetry);
      if (!assign(line, ob_1)) return n;

      const FT delta = CGAL::sqrt(FT(1) - cos_symmetry * cos_symmetry);
      const Point projected_origin = line.projection(CGAL::ORIGIN);
      Vector line_vector(line);
      line_vector /= CGAL::sqrt(line_vector * line_vector);
      const Point pt1 = projected_origin + delta * line_vector;
      const Point pt2 = projected_origin - delta * line_vector;

      if (CGAL::squared_distance(pt_normal, pt1) <= CGAL::squared_distance(pt_normal, pt2))
        return Vector(CGAL::ORIGIN, pt1);
      else
        return Vector(CGAL::ORIGIN, pt2);
    } else return n;
  }

  template<typename Traits>
  typename Traits::Vector_3 regularize_normals_from_prior(
    const typename Traits::Vector_3& np,
    const typename Traits::Vector_3& n,
    const typename Traits::Vector_3& symmetry_direction,
    const typename Traits::FT cos_symmetry) {

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Line_3 Line;
    typedef typename Traits::Plane_3 Plane;

    const Plane plane_orthogonality(CGAL::ORIGIN, np);
    const Point pt_symmetry = CGAL::ORIGIN + cos_symmetry * symmetry_direction;
    const Plane plane_symmetry(pt_symmetry, symmetry_direction);

    Line line;
    const CGAL::Object ob_1 = CGAL::intersection(plane_orthogonality, plane_symmetry);
    if (!assign(line, ob_1))
      return regularize_normal<Traits>(n, symmetry_direction, cos_symmetry);

    const Point projected_origin = line.projection(CGAL::ORIGIN);
    const FT R = CGAL::squared_distance(Point(CGAL::ORIGIN), projected_origin);

    if (R <= 1) { // 2 (or 1) possible points intersecting the unit sphere and line
      const FT delta = std::sqrt (FT(1) - R);
      Vector line_vector(line);
      line_vector /= CGAL::sqrt(line_vector * line_vector);
      const Point pt1 = projected_origin + delta * line_vector;
      const Point pt2 = projected_origin - delta * line_vector;

      const Point pt_n = CGAL::ORIGIN + n;
      if (CGAL::squared_distance(pt_n, pt1) <= CGAL::squared_distance(pt_n, pt2))
        return Vector(CGAL::ORIGIN, pt1);
      else
        return Vector(CGAL::ORIGIN, pt2);
    } else // no point intersecting the unit sphere and line
      return regularize_normal<Traits>(n, symmetry_direction, cos_symmetry);
  }

  template<
  typename Traits,
  typename PointRange,
  typename PointMap,
  typename IndexMap>
  void compute_centroids_and_areas(
    const PointRange& points,
    PointMap point_map,
    const std::size_t nb_planes,
    IndexMap index_map,
    std::vector<typename Traits::Point_3>& centroids,
    std::vector<typename Traits::FT>& areas) {

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;

    std::vector< std::vector<Point> > listp(nb_planes);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const int idx = get(index_map, i);
      if (idx != -1)
        listp[std::size_t(idx)].push_back(
          get(point_map, *(points.begin() + i)));
    }

    centroids.reserve(nb_planes);
    areas.reserve(nb_planes);
    for (std::size_t i = 0; i < nb_planes; ++i) {
      centroids.push_back(
        CGAL::centroid(listp[i].begin(), listp[i].end()));
      areas.push_back(FT(listp[i].size() / FT(100)));
    }
  }

  template<
  typename Traits,
  typename PlaneRange,
  typename PlaneMap>
  void compute_parallel_clusters(
    const PlaneRange& planes,
    PlaneMap plane_map,
    std::vector<Plane_cluster<Traits> >& clusters,
    const std::vector<typename Traits::FT>& areas,
    const typename Traits::FT tolerance_cosangle,
    const typename Traits::Vector_3& symmetry_direction) {

    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector;

    // Find pairs of epsilon-parallel primitives and store them in parallel_planes.
    std::vector< std::vector<std::size_t> > parallel_planes(planes.size());
    for (std::size_t i = 0; i < std::size_t(planes.size()); ++i) {
      const auto it = planes.begin() + i;
      const Vector v1 = get(plane_map, *it).orthogonal_vector();

      for (std::size_t j = 0; j < std::size_t(planes.size()); ++j) {
        if (i == j) continue;

        const auto it2 = planes.begin() + j;
        const Vector v2 = get(plane_map, *it2).orthogonal_vector();

        if (CGAL::abs(v1 * v2) > FT(1) - tolerance_cosangle)
          parallel_planes[i].push_back(j);
      }
    }

    std::vector<bool> is_available(planes.size(), true);
    for (std::size_t i = 0; i < std::size_t(planes.size()); ++i) {
      if (is_available[i]) {
        const auto& plane = get(plane_map, *(planes.begin() + i));
        is_available[i] = false;

        clusters.push_back(Plane_cluster<Traits>());
        Plane_cluster<Traits>& clu = clusters.back();

        // Initialize containers.
        clu.planes.push_back(i);

        std::vector<std::size_t> index_container_former_ring_parallel;
        index_container_former_ring_parallel.push_back(i);
        std::list<std::size_t> index_container_current_ring_parallel;

        // Propagate over the pairs of epsilon-parallel primitives.
        bool propagation = true;
        clu.normal = plane.orthogonal_vector();
        clu.area = areas[i];
        do {
          propagation = false;

          for (std::size_t k = 0; k < index_container_former_ring_parallel.size(); ++k) {
            const std::size_t plane_index = index_container_former_ring_parallel[k];
            for (std::size_t l = 0; l < parallel_planes[plane_index].size(); ++l) {
              const std::size_t it = parallel_planes[plane_index][l];

              Vector normal_it =
                get(plane_map, *(planes.begin() + it)).orthogonal_vector();
              if (is_available[it] &&
                CGAL::abs(normal_it * clu.normal) > FT(1) - tolerance_cosangle) {

                propagation = true;
                index_container_current_ring_parallel.push_back(it);
                is_available[it] = false;

                if (clu.normal * normal_it < FT(0))
                  normal_it = -normal_it;

                clu.normal = FT(clu.area) * clu.normal + FT(areas[it]) * normal_it;
                const FT norm = FT(1) / CGAL::sqrt(clu.normal.squared_length());
                clu.normal = norm * clu.normal;
                clu.area += areas[it];
              }
            }
          }

          // Update containers.
          index_container_former_ring_parallel.clear();
          for (auto it = index_container_current_ring_parallel.begin();
          it != index_container_current_ring_parallel.end(); ++it) {
            index_container_former_ring_parallel.push_back(*it);
            clu.planes.push_back(*it);
          }
          index_container_current_ring_parallel.clear();

        } while (propagation);

        if (symmetry_direction != CGAL::NULL_VECTOR) {
          clu.cosangle_symmetry = symmetry_direction * clu.normal;
          if (clu.cosangle_symmetry < FT(0)) {
            clu.normal = -clu.normal;
            clu.cosangle_symmetry = -clu.cosangle_symmetry;
          }
        }
      }
    }
    is_available.clear();
  }

  template<typename Traits>
  void cluster_symmetric_cosangles(
    std::vector<Plane_cluster<Traits> >& clusters,
    const typename Traits::FT tolerance_cosangle,
    const typename Traits::FT tolerance_cosangle_ortho) {

    typedef typename Traits::FT FT;

    std::vector<FT> cosangle_centroids;
    std::vector<std::size_t> list_cluster_index;
    for (std::size_t i = 0; i < clusters.size(); ++i)
      list_cluster_index.push_back(static_cast<std::size_t>(-1));

    std::size_t mean_index = 0;
    for (std::size_t i = 0; i < clusters.size(); ++i) {
      if (list_cluster_index[i] == static_cast<std::size_t>(-1)) {
        list_cluster_index[i] = mean_index;
        FT mean = clusters[i].area * clusters[i].cosangle_symmetry;
        FT mean_area = clusters[i].area;

        for (std::size_t j = i + 1; j < clusters.size(); ++j) {
          if (list_cluster_index[j] == static_cast<std::size_t>(-1) &&
            CGAL::abs(clusters[j].cosangle_symmetry - mean / mean_area) < tolerance_cosangle_ortho) {

            list_cluster_index[j] = mean_index;
            mean_area += clusters[j].area;
            mean += clusters[j].area * clusters[j].cosangle_symmetry;
          }
        }
        ++mean_index;
        mean /= mean_area;
        cosangle_centroids.push_back(mean);
      }
    }

    for (std::size_t i = 0; i < cosangle_centroids.size(); ++i) {
      if (cosangle_centroids[i] < tolerance_cosangle_ortho)
        cosangle_centroids[i] = FT(0);
      else if (cosangle_centroids[i] > FT(1) - tolerance_cosangle)
        cosangle_centroids[i] = FT(1);
    }

    for (std::size_t i = 0; i < clusters.size(); ++i)
      clusters[i].cosangle_symmetry = cosangle_centroids[list_cluster_index[i]];
  }

  template<typename Traits>
  void subgraph_mutually_orthogonal_clusters(
    std::vector< Plane_cluster<Traits> >& clusters,
    const typename Traits::Vector_3& symmetry_direction) {

    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector;

    std::vector< std::vector<std::size_t> > subgraph_clusters;
    std::vector<std::size_t> subgraph_clusters_max_area_index;

    for (std::size_t i = 0; i < clusters.size(); ++i)
      clusters[i].is_free = true;

    for (std::size_t i = 0; i < clusters.size(); ++i) {
      if (clusters[i].is_free) {
        clusters[i].is_free = false;
        FT max_area = clusters[i].area;
        std::size_t index_max_area = i;

        // Initialize containers.
        std::vector<std::size_t> index_container;
        index_container.push_back(i);
        std::vector<std::size_t> index_container_former_ring;
        index_container_former_ring.push_back(i);
        std::list<std::size_t> index_container_current_ring;

        // Propagate.
        bool propagation = true;
        do {
          propagation = false;

          // Neighbors.
          for (std::size_t k = 0; k < index_container_former_ring.size(); ++k) {
            const std::size_t cluster_index_1 = index_container_former_ring[k];

            for (std::size_t j = 0; j < clusters[cluster_index_1].orthogonal_clusters.size(); ++j) {
              const std::size_t cluster_index_2 = clusters[cluster_index_1].orthogonal_clusters[j];
              if (clusters[cluster_index_2].is_free) {
                propagation = true;
                index_container_current_ring.push_back(cluster_index_2);
                clusters[cluster_index_2].is_free = false;

                if (max_area < clusters[cluster_index_2].area) {
                  max_area = clusters[cluster_index_2].area;
                  index_max_area = cluster_index_2;
                }
              }
            }
          }

          // Update containers.
          index_container_former_ring.clear();
          for (auto it = index_container_current_ring.begin();
            it != index_container_current_ring.end(); ++it) {

            index_container_former_ring.push_back(*it);
            index_container.push_back(*it);
          }
          index_container_current_ring.clear();

        } while (propagation);

        subgraph_clusters.push_back(index_container);
        subgraph_clusters_max_area_index.push_back(index_max_area);
      }
    }

    // Create subgraphs of mutually orthogonal clusters in which the
    // largest cluster is excluded and then store them in subgraph_clusters_prop.
    std::vector< std::vector<std::size_t> > subgraph_clusters_prop;
    for (std::size_t i = 0; i < subgraph_clusters.size(); ++i) {

      const std::size_t index = subgraph_clusters_max_area_index[i];
      std::vector<std::size_t> subgraph_clusters_prop_temp;
      for (std::size_t j = 0; j < subgraph_clusters[i].size(); ++j)
        if (subgraph_clusters[i][j] != index)
          subgraph_clusters_prop_temp.push_back(subgraph_clusters[i][j]);
      subgraph_clusters_prop.push_back(subgraph_clusters_prop_temp);
    }

    // Regularize cluster normals : in each subgraph, we start
    // from the largest area cluster and we propagate over the subgraph
    // by regularizing the normals of the clusters according to the
    // orthogonality and cos angle to symmetry direction.
    for (std::size_t i = 0; i < clusters.size(); ++i)
      clusters[i].is_free = true;

    for (std::size_t i = 0; i < subgraph_clusters_prop.size(); ++i) {
      const std::size_t index_current = subgraph_clusters_max_area_index[i];
      const Vector vec_current = regularize_normal<Traits>(
        clusters[index_current].normal,
        symmetry_direction,
        clusters[index_current].cosangle_symmetry);

      clusters[index_current].normal = vec_current;
      clusters[index_current].is_free = false;

      // Initialize containers.
      std::vector<std::size_t> index_container;
      index_container.push_back(index_current);
      std::vector<std::size_t> index_container_former_ring;
      index_container_former_ring.push_back(index_current);
      std::list<std::size_t> index_container_current_ring;

      // Propagate.
      bool propagation = true;
      do {
        propagation = false;

        // Neighbors.
        for (std::size_t k = 0; k < index_container_former_ring.size(); ++k) {
          const std::size_t cluster_index_1 = index_container_former_ring[k];

          for (std::size_t j = 0; j < clusters[cluster_index_1].orthogonal_clusters.size(); ++j) {
            const std::size_t cluster_index_2 = clusters[cluster_index_1].orthogonal_clusters[j];
            if (clusters[cluster_index_2].is_free) {
              propagation = true;
              index_container_current_ring.push_back(cluster_index_2);
              clusters[cluster_index_2].is_free = false;

              const Vector new_vect = regularize_normals_from_prior<Traits>(
                clusters[cluster_index_1].normal,
                clusters[cluster_index_2].normal,
                symmetry_direction,
                clusters[cluster_index_2].cosangle_symmetry);

              clusters[cluster_index_2].normal = new_vect;
            }
          }
        }

        // Update containers.
        index_container_former_ring.clear();
        for (auto it = index_container_current_ring.begin();
          it != index_container_current_ring.end(); ++it) {

          index_container_former_ring.push_back(*it);
          index_container.push_back(*it);
        }
        index_container_current_ring.clear();
      } while (propagation);
    }
  }

} // internal
} // Shape_regularization
} // CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_UTILS_H
