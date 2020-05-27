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

#ifndef CGAL_SHAPE_REGULARIZATION_REGULARIZE_SEGMENTS_H
#define CGAL_SHAPE_REGULARIZATION_REGULARIZE_SEGMENTS_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/QP_regularization.h>
#include <CGAL/Shape_regularization/Solvers/OSQP_quadratic_program.h>

#include <CGAL/Shape_regularization/Segments/Angle_regularization_2.h>
#include <CGAL/Shape_regularization/Segments/Offset_regularization_2.h>
#include <CGAL/Shape_regularization/Segments/Delaunay_neighbor_query_2.h>

#include <CGAL/Shape_regularization/Segments/Parallel_groups_2.h>
#include <CGAL/Shape_regularization/Segments/Orthogonal_groups_2.h>
#include <CGAL/Shape_regularization/Segments/Collinear_groups_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  template<
  typename InputRange,
  typename NeighborQuery,
  typename RegularizationType,
  typename QPSolver,
  typename GeomTraits>
  void regularize_segments(
    const InputRange& input_range,
    NeighborQuery& neighbor_query,
    RegularizationType& regularization_type,
    QPSolver& quadratic_program,
    GeomTraits traits) {

    using Regularizer =
      CGAL::Shape_regularization::QP_regularization<
      GeomTraits, InputRange, NeighborQuery, RegularizationType, QPSolver>;

    Regularizer regularizer(
      input_range, neighbor_query, regularization_type, quadratic_program, traits);
    regularizer.regularize();
  }

  template<
  typename InputRange,
  typename NeighborQuery,
  typename RegularizationType,
  typename GeomTraits>
  void regularize_segments(
    const InputRange& input_range,
    NeighborQuery& neighbor_query,
    RegularizationType& regularization_type,
    GeomTraits traits) {

    using FT = typename GeomTraits::FT;
    using QP = CGAL::Shape_regularization::OSQP_quadratic_program<FT>;
    QP quadratic_program;

    regularize_segments(
      input_range, neighbor_query, regularization_type, quadratic_program, traits);
  }

  template<
  typename InputRange,
  typename NeighborQuery,
  typename RegularizationType>
  void regularize_segments(
    const InputRange& input_range,
    NeighborQuery& neighbor_query,
    RegularizationType& regularization_type) {

    using Iterator_type = typename InputRange::const_iterator;
    using Segment_2 = typename std::iterator_traits<Iterator_type>::value_type;
    using GeomTraits = typename Kernel_traits<Segment_2>::Kernel;
    GeomTraits traits;

    using FT = typename GeomTraits::FT;
    using QP = CGAL::Shape_regularization::OSQP_quadratic_program<FT>;
    QP quadratic_program;

    regularize_segments(
      input_range, neighbor_query, regularization_type, quadratic_program, traits);
  }

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_REGULARIZE_SEGMENTS_H
