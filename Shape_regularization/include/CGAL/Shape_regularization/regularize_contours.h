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
// Author(s)     : Dmitry Anisimov, Simon Giraudot
//

#ifndef CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H
#define CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/internal/Contour_regularization_2.h>

#include <CGAL/Shape_regularization/Contours/Longest_direction_2.h>
#include <CGAL/Shape_regularization/Contours/Multiple_directions_2.h>
#include <CGAL/Shape_regularization/Contours/User_defined_directions_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes closed contours.

    Given a set of counterclockwise ordered 2D points connected by segments, which form a closed contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam ContourDirections
    must be a model of `ContourDirections`.

    \tparam OutputIterator
    must be an output iterator whose value type is `GeomTraits::Point_2`.

    \tparam NamedParameters
    a sequence of \ref sr_namedparameters "Named Parameters".

    \tparam PointMap
    must be a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    optional sequence of \ref sr_namedparameters "Named Parameters"
    among the ones listed below

    \param point_map
    an instance of `PointMap`

    \param traits
    an instance of `GeomTraits`

    \cgalNamedParamsBegin
      \cgalParamBegin{max_offset}
        max distance between two parallel and consecutive contour edges,
        the default is 0.5 unit length
      \cgalParamEnd
    \cgalNamedParamsEnd

    \return an output iterator.

    \pre input_range.size() >= 3
  */
  template<
  typename InputRange,
  typename ContourDirections,
  typename OutputIterator,
  typename NamedParameters,
  typename PointMap,
  typename GeomTraits>
  OutputIterator regularize_closed_contour(
    const InputRange& input_range,
    const ContourDirections& directions,
    OutputIterator contour,
    const NamedParameters np,
    const PointMap point_map,
    const GeomTraits traits) {

    CGAL_precondition(input_range.size() >= 3);
    using Contour_regularizer =
    CGAL::Shape_regularization::internal::Contour_regularization_2<
      CGAL::Shape_regularization::internal::CLOSED,
      ContourDirections, GeomTraits>;

    Contour_regularizer regularizer(
      directions, input_range, point_map, np, traits);
    return regularizer.regularize(contour);
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes closed contours.

    Given a set of counterclockwise ordered 2D points connected by segments, which form a closed contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    This function infers a traits class `GeomTraits` from the `InputRange` iterator's value type.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam ContourDirections
    must be a model of `ContourDirections`.

    \tparam OutputIterator
    must be an output iterator whose value type is `GeomTraits::Point_2`.

    \tparam NamedParameters
    a sequence of \ref sr_namedparameters "Named Parameters".

    \tparam PointMap
    must be a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`. %Default is the
    `CGAL::Identity_property_map`.

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    optional sequence of \ref sr_namedparameters "Named Parameters"
    among the ones listed below

    \param point_map
    an instance of `PointMap`, if not provided, the default is used

    \cgalNamedParamsBegin
      \cgalParamBegin{max_offset}
        max distance between two parallel and consecutive contour edges,
        the default is 0.5 unit length
      \cgalParamEnd
    \cgalNamedParamsEnd

    \return an output iterator.

    \pre input_range.size() >= 3
  */
  template<
  typename InputRange,
  typename ContourDirections,
  typename OutputIterator,
  typename NamedParameters,
  typename PointMap = CGAL::Identity_property_map<
  typename std::iterator_traits< typename InputRange::const_iterator >::value_type > >
  OutputIterator regularize_closed_contour(
    const InputRange& input_range,
    const ContourDirections& directions,
    OutputIterator contour,
    const NamedParameters np,
    const PointMap point_map = PointMap()) {

    CGAL_precondition(input_range.size() >= 3);
    using Iterator_type = typename InputRange::const_iterator;
    using Point_2 = typename std::iterator_traits<Iterator_type>::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    GeomTraits traits;

    return regularize_closed_contour(
      input_range, directions, contour, np, point_map, traits);
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes open contours.

    Given a set of counterclockwise ordered 2D points connected by segments, which form an open contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam ContourDirections
    must be a model of `ContourDirections`.

    \tparam OutputIterator
    must be an output iterator whose value type is `GeomTraits::Point_2`.

    \tparam NamedParameters
    a sequence of \ref sr_namedparameters "Named Parameters".

    \tparam PointMap
    must be a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    optional sequence of \ref sr_namedparameters "Named Parameters"
    among the ones listed below

    \param point_map
    an instance of `PointMap`

    \param traits
    an instance of `GeomTraits`

    \cgalNamedParamsBegin
      \cgalParamBegin{max_offset}
        max distance between two parallel and consecutive contour edges,
        the default is 0.5 unit length
      \cgalParamEnd
    \cgalNamedParamsEnd

    \return an output iterator.

    \pre input_range.size() >= 2
  */
  template<
  typename InputRange,
  typename ContourDirections,
  typename OutputIterator,
  typename NamedParameters,
  typename PointMap,
  typename GeomTraits>
  OutputIterator regularize_open_contour(
    const InputRange& input_range,
    const ContourDirections& directions,
    OutputIterator contour,
    const NamedParameters np,
    const PointMap point_map,
    const GeomTraits traits) {

    CGAL_precondition(input_range.size() >= 2);
    using Contour_regularizer =
    CGAL::Shape_regularization::internal::Contour_regularization_2<
      CGAL::Shape_regularization::internal::OPEN,
      ContourDirections, GeomTraits>;

    Contour_regularizer regularizer(
      directions, input_range, point_map, np, traits);
    return regularizer.regularize(contour);
  }

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief regularizes open contours.

    Given a set of counterclockwise ordered 2D points connected by segments, which form an open contour,
    this function enables to reinforce three types of regularities among consecutive edges of this contour:
    - *Parallelism*: contour edges, which are detected as near parallel, are made exactly parallel.
    - *Orthogonality*: contour edges, which are detected as near orthogonal, are made exactly orthogonal.
    - *Collinearity*: parallel contour edges, which are detected as near collinear, are made exactly collinear.

    The principal directions of the contour are provided via the concept `ContourDirections`.

    This function infers a traits class `GeomTraits` from the `InputRange` iterator's value type.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam ContourDirections
    must be a model of `ContourDirections`.

    \tparam OutputIterator
    must be an output iterator whose value type is `GeomTraits::Point_2`.

    \tparam NamedParameters
    a sequence of \ref sr_namedparameters "Named Parameters".

    \tparam PointMap
    must be a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`. %Default is the
    `CGAL::Identity_property_map`.

    \param input_range
    a const range of ordered points, which form a contour

    \param directions
    estimated contour directions towards which the contour edges are oriented

    \param contour
    an output iterator with points of the regularized contour

    \param np
    optional sequence of \ref sr_namedparameters "Named Parameters"
    among the ones listed below

    \param point_map
    an instance of `PointMap`, if not provided, the default is used

    \cgalNamedParamsBegin
      \cgalParamBegin{max_offset}
        max distance between two parallel and consecutive contour edges,
        the default is 0.5 unit length
      \cgalParamEnd
    \cgalNamedParamsEnd

    \return an output iterator.

    \pre input_range.size() >= 2
  */
  template<
  typename InputRange,
  typename ContourDirections,
  typename OutputIterator,
  typename NamedParameters,
  typename PointMap = CGAL::Identity_property_map<
  typename std::iterator_traits< typename InputRange::const_iterator >::value_type > >
  OutputIterator regularize_open_contour(
    const InputRange& input_range,
    const ContourDirections& directions,
    OutputIterator contour,
    const NamedParameters np,
    const PointMap point_map = PointMap()) {

    CGAL_precondition(input_range.size() >= 2);
    using Iterator_type = typename InputRange::const_iterator;
    using Point_2 = typename std::iterator_traits<Iterator_type>::value_type;
    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    GeomTraits traits;

    return regularize_open_contour(
      input_range, directions, contour, np, point_map, traits);
  }

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_REGULARIZE_CONTOURS_H
