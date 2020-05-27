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

#ifndef CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H
#define CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Contour_base_2.h>

namespace CGAL {
namespace Shape_regularization {
namespace Contours {

  /*!
    \ingroup PkgShapeRegularizationRefContours

    \brief Estimates possibly multiple principal directions of the contour
    based on the user-defined min length and max angle bounds.

    This algorithm finds the best-fit edges of the contour with respect to the
    user-defined conditions and sets them as its principal directions.

    \tparam GeomTraits
    must be a model of `Kernel`.

    \tparam InputRange
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam PointMap
    must be a model of `ReadablePropertyMap` whose key type is the value type of the input
    range and value type is `GeomTraits::Point_2`. %Default is the
    `CGAL::Identity_property_map<typename GeomTraits::Point_2>`.

    \cgalModels `ContourDirections`
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Multiple_directions_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Segment_2 = typename Traits::Segment_2;
    using Direction_2 = typename Traits::Direction_2;

    using FT_pair = std::pair<FT, FT>;
    using Base = internal::Contour_base_2<Traits>;
    using Segment_wrapper_2 = typename Base::Segment_wrapper_2;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref sr_namedparameters "Named Parameters".

      \param input_range
      a const range of points, which form a contour

      \param np
      optional sequence of \ref sr_namedparameters "Named Parameters"
      among the ones listed below

      \param is_closed
      indicates whether the contour is closed or open

      \param point_map
      an instance of `PointMap`, if not provided, the default is used

      \cgalNamedParamsBegin
        \cgalParamBegin{max_angle}
          max angle deviation in degrees between a contour edge and a given
          principal direction, the default is 25 degrees
        \cgalParamEnd
        \cgalParamBegin{min_length}
          min length in meters of a contour edge whose direction can be taken
          as a principal direction, the default is 3 meters
        \cgalParamEnd
      \cgalNamedParamsEnd

      \pre `input_range.size() >= 3` for closed contours
      \pre `input_range.size() >= 2` for open contours
    */
    template<typename NamedParameters>
    Multiple_directions_2(
      const InputRange& input_range,
      const bool is_closed,
      const NamedParameters np,
      const PointMap point_map = PointMap()) :
    m_input_range(input_range),
    m_point_map(point_map) {

      CGAL_precondition(input_range.size() >= 2);
      m_min_length_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_length), FT(3));
      m_max_angle_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(25));

      CGAL_precondition(
        m_min_length_2 >= FT(0));
      CGAL_precondition(
        m_max_angle_2 >= FT(0) &&
        m_max_angle_2 <= FT(90));

      if (is_closed)
        estimate_closed(m_bounds, m_directions, m_assigned);
      else
        estimate_open(m_bounds, m_directions, m_assigned);

      if (verbose()) {
        std::cout << "* assigned directions: ";
        for (std::size_t direction_index : m_assigned)
          std::cout << direction_index << " ";
        std::cout << std::endl;
      }
    }

    /// @}

    /// \name Directions
    /// @{

    /*!
      \brief orients a given `segment` with the index `query_index` with respect
      to the best found principal direction.

      \param query_index an index of the `segment` in the input contour, in other words,
      the segment's source point is the point in the contour with the index `query_index`

      \param segment a segment to be oriented

      \pre `query_index >= 0 && query_index < input_range.size()` for closed contours
      \pre `query_index >= 0 && query_index < input_range.size() - 1` for open contours
    */
    void orient(
      const std::size_t query_index,
      Segment_2& segment) const {

      m_base.apply_rotation_to_segment(
        m_bounds, m_directions, m_assigned,
        query_index, segment);
    }

    /// @}

    /// \name Miscellaneous
    /// @{

    /*!
      \brief returns the number of principal directions of the contour.

      The returned number is always greater or equal than one.
    */
    const std::size_t number_of_directions() const {
      return m_directions.size();
    }

    /// @}

    // EXTRA METHODS TO TEST THE CLASS!
    /// \cond SKIP_IN_MANUAL
    const std::vector<Direction_2>& get_directions() const {
      return m_directions;
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Point_map m_point_map;
    const Base m_base;

    FT m_min_length_2;
    FT m_max_angle_2;

    std::vector<FT_pair> m_bounds;
    std::vector<Direction_2> m_directions;
    std::vector<std::size_t> m_assigned;

    const bool verbose() const {
      return m_base.verbose();
    }

    void estimate_closed(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      std::vector<Segment_wrapper_2> wraps;
      m_base.initialize_closed(
        m_min_length_2, m_input_range, m_point_map, wraps);
      m_base.estimate_initial_directions(
        m_max_angle_2, wraps, bounds, directions, assigned);

      if (directions.size() <= 1) {
        m_base.set_longest_direction(
          wraps, bounds, directions, assigned);
      } else {
        m_base.unify_along_contours_closed(wraps, assigned);
        m_base.correct_directions_closed(wraps, assigned);
        m_base.readjust_directions(wraps, assigned, directions);
      }
    }

    void estimate_open(
      std::vector<FT_pair>& bounds,
      std::vector<Direction_2>& directions,
      std::vector<std::size_t>& assigned) const {

      std::vector<Segment_wrapper_2> wraps;
      m_base.initialize_open(
        m_min_length_2, m_input_range, m_point_map, wraps);
      m_base.estimate_initial_directions(
        m_max_angle_2, wraps, bounds, directions, assigned);

      if (directions.size() <= 1) {
        m_base.set_longest_direction(
          wraps, bounds, directions, assigned);
      } else {
        m_base.unify_along_contours_open(wraps, assigned);
        m_base.correct_directions_open(wraps, assigned);
        m_base.readjust_directions(wraps, assigned, directions);
      }
    }
  };

} // namespace Contours
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_MULTIPLE_PRINCIPAL_DIRECTIONS_2_H
