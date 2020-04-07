// Copyright (c) 2020 GeometryFactory Sarl (France).
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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Andreas Fabri
//

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <vector>
#include <utility>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

// Internal includes.
#include <CGAL/Shape_regularization/enum.h>
#include <CGAL/Shape_regularization/internal/Closed_contour_regularization_2.h>
#include <CGAL/Shape_regularization/internal/Open_contour_regularization_2.h>

// TODO:
// * What about parameterizing this class by the Direction_estimator class: we have one 
// for the longest directions, one for the multiple, and one for the user-defined.
// All these estimators should return m_directions, m_assigned, and m_bounds.

namespace CGAL {
namespace Shape_regularization {

  struct OPEN { };
  struct CLOSED { };

  namespace internal {
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_iterator, iterator, false)
  }

  template<
  typename InputRange, 
  typename NamedParameters,
  bool has_nested_iterator=internal::Has_nested_type_iterator<InputRange>::value>
  class GetPointMap {
    
    typedef typename std::iterator_traits<typename InputRange::iterator>::value_type Point;
    typedef typename CGAL::Identity_property_map<Point> DefaultPMap;

  public:
    typedef typename internal_np::Lookup_named_param_def<
    internal_np::point_t,
    NamedParameters,
    DefaultPMap
    > ::type  type;

    typedef typename internal_np::Lookup_named_param_def<
    internal_np::point_t,
    NamedParameters,
    DefaultPMap
    > ::type  const_type;
  };

  // To please compiler instantiating non valid overloads.
  template<
  typename InputRange, 
  typename NamedParameters>
  class GetPointMap<InputRange, NamedParameters, false> {
    struct Dummy_point { };
  public:
    typedef typename CGAL::Identity_property_map<Dummy_point> type;
    typedef typename CGAL::Identity_property_map<Dummy_point> const_type;
  };

  /*!
    \ingroup PkgShapeRegularizationRef
    
    \brief Contour regularization algorithm.

    This algorithm enables to regularize both open and closed contours.

    \tparam GeomTraits 
    must be a model of `Kernel`.
    
    \tparam InputRange 
    must be a model of `ConstRange`.

    \tparam PointMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Point_2`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename Tag>
  class Contour_regularization_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    
    using FT = typename Traits::FT;

    using Regularization = typename std::conditional<
      std::is_same<Tag, CLOSED>::value,
      internal::Closed_contour_regularization_2<Traits>,
      internal::Open_contour_regularization_2<Traits> >::type;

    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam
      a sequence of \ref pmp_namedparameters "Named Parameters"

      \param input_range
      a range of points, which form a contour

      \param np
      optional sequence of \ref pmp_namedparameters "Named Parameters" 
      among the ones listed below

      \pre `input_range.size() >= 3` for closed contours
      \pre `input_range.size() >= 2` for open contours
    */
    template<typename NamedParameters>
    Contour_regularization_2(
      Input_range& input_range,
      const NamedParameters& np) { 

      using Point_map = typename GetPointMap<InputRange, NamedParameters>::type;
      const Point_map point_map = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::point_map), Point_map());
      const FT min_length_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::min_length), FT(3));
      const FT max_angle_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(25));
      const FT max_offset_2 = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_offset), FT(1) / FT(2));

      m_regularization = std::make_shared<Regularization>(
        min_length_2, max_angle_2, max_offset_2);
      m_regularization->initialize(input_range, point_map);

      if (m_regularization->verbose()) {
        std::cout << "* parameters: " << std::endl;
        std::cout << "- min length 2: " << min_length_2 << std::endl;
        std::cout << "- max angle 2: " << max_angle_2 << std::endl;
        std::cout << "- max offset 2: " << max_offset_2 << std::endl;
      }
    }

    /// @}

    /// \name Functions
    /// @{

    /*!
      \brief sets principal directions of the contour.

      This method sets the user-defined principal contour directions. All contour 
      edges will be regularized with respect to these directions.

      \tparam DirectionRange
      must be a model of `ConstRange`.

      \tparam DirectionMap
      must be an `LvaluePropertyMap` whose key type is the value type of the `DirectionRange`
      and value type is `GeomTraits::Direction_2`.

      \param direction_range 
      a range of user-defined directions

      \param direction_map
      an instance of `DirectionMap`
    */
    template<
    typename DirectionRange,
    typename DirectionMap>
    void set_principal_directions(
      const DirectionRange& direction_range,
      const DirectionMap direction_map) {

      m_regularization->set_principal_directions(
        direction_range, direction_map);
    }

    /*!
      \brief estimates principal directions of the contour.

      This method estimates the principal contour directions automatically.

      \param direction_type
      indicates which type of principal directions should be used, 
      the default is `Direction_type::LONGEST`
    */
    void estimate_principal_directions(
      const Direction_type direction_type = Direction_type::LONGEST) {
      m_regularization->estimate_principal_directions(
        direction_type);
    }

    /*!
      \brief executes the contour regularization algorithm.

      This method regularizes the contour with respect to the defined principal directions.
      That is it sets all other not principal segments either orthogonal or collinear 
      to the chosen directions.

      \tparam OutputIterator 
      must be an output iterator whose value type is `GeomTraits::Point_2`.

      \param contour
      an `OutputIterator` with contour points
    */
    template<typename OutputIterator>
    void regularize(
      OutputIterator contour) {
      m_regularization->regularize(
        contour);
    }

    /*!
      \brief returns the number of principal directions in the contour.
    */
    const std::size_t number_of_principal_directions() const {
      return m_regularization->number_of_principal_directions();
    }

  private:
    std::shared_ptr<Regularization> m_regularization;
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
