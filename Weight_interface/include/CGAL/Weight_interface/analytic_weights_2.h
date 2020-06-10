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

#ifndef CGAL_GENERALIZED_ANALYTIC_WEIGHTS_2_H
#define CGAL_GENERALIZED_ANALYTIC_WEIGHTS_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/Generalized_weights_2/enum_2.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Wachspress_weights_2.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Discrete_harmonic_weights_2.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Mean_value_weights_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `AnalyticTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Wachspress_weights_2<PointRange, GeomTraits> wachspress(
      polygon, policy, traits);
    return wachspress(query, weights);
  }

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D Wachspress weights.

    This function computes 2D Wachspress weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Wachspress_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator wachspress_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator weights,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return wachspress_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `AnalyticTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Discrete_harmonic_weights_2<PointRange, GeomTraits> discrete_harmonic(
      polygon, policy, traits);
    return discrete_harmonic(query, weights);
  }

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D discrete harmonic weights.

    This function computes 2D discrete harmonic weights at a given `query` point
    with respect to the vertices of a strictly convex `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Discrete_harmonic_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a strictly convex polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
    \pre `polygon is strictly convex`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator discrete_harmonic_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator weights,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return discrete_harmonic_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `GeomTraits::FT`.

    \tparam GeomTraits
    is a model of `AnalyticTraits_2`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param traits
    An instance of `GeomTraits`.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename OutputIterator,
  typename GeomTraits>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const typename GeomTraits::Point_2& query,
    OutputIterator weights,
    const GeomTraits traits,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    Mean_value_weights_2<PointRange, GeomTraits> mean_value(
      polygon, policy, traits);
    return mean_value(query, weights);
  }

  /*!
    \ingroup PkgWeightInterfaceRef2DFunctions

    \brief computes 2D mean value weights.

    This function computes 2D mean value weights at a given `query` point
    with respect to the vertices of a simple `polygon`, that is one
    weight per vertex. The weights are returned in `weights`.

    Internally, the class `CGAL::Generalized_weights::Mean_value_weights_2` is used.
    If one needs a flexible API, please refer to that class. If you want to handle
    multiple query points, you better use that class, too. When using this function,
    internal memory is allocated for each query point, while when using the class,
    it is allocated only once, which is much more efficient.

    This function infers a traits class from the `Point_2` class.

    \tparam PointRange
    is a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam OutputIterator
    is an output iterator whose value type is `Kernel_traits<Point_2>::Kernel::FT`.

    \param polygon
    An instance of `PointRange` with 2D points, which form a simple polygon.

    \param query
    A query point.

    \param weights
    An output iterator that stores the computed weights.

    \param policy
    One of the `CGAL::Generalized_weights::Computation_policy_2`.
    The default is `CGAL::Generalized_weights::Computation_policy_2::DEFAULT`.

    \return an output iterator.

    \pre `polygon.size() >= 3`
    \pre `polygon is simple`
  */
  template<
  typename PointRange,
  typename Point_2,
  typename OutputIterator>
  OutputIterator mean_value_weights_2(
    const PointRange& polygon,
    const Point_2& query,
    OutputIterator weights,
    const Computation_policy_2 policy =
    Computation_policy_2::DEFAULT) {

    using GeomTraits = typename Kernel_traits<Point_2>::Kernel;
    return mean_value_weights_2(
      polygon, query, weights, GeomTraits(), policy);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_ANALYTIC_WEIGHTS_2_H
