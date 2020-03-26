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
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Andreas Fabri
//

#ifndef CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/number_utils.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef
    
    \brief Contour regularization algorithm.

    This algorithm enables to regularize both open and closed contours.

    \tparam GeomTraits 
    must be a model of `Kernel`.
    
    \tparam InputRange 
    must be a model of `ConstRange`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Segment_2`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename SegmentMap>
  class Contour_regularization_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param parameters
      input parameters
    */
    Contour_regularization_2(
      std::string parameters) :
    m_parameters(parameters) { 
      
    }

    /// @}

    /// \name Functions
    /// @{

    /*!
      \brief sets principal directions of the contour.

      This method sets the user-defined principal contour directions.

      \param input_range
      a range of input segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` 
      to `GeomTraits::Segment_2`

      \param directions 
      a set of indices where each index represents a segment from the `input_range` 
      that is chosen as one of the principal directions

      \pre `input_range.size() > 1`
      \pre `directions.size() == input_range.size()`
    */
    void set_principal_directions(
      const InputRange& input_range,
      const SegmentMap segment_map,
      const std::vector<std::size_t>& directions) {

      CGAL_precondition(input_range.size() > 1);
      CGAL_precondition(directions.size() == input_range.size());
    }

    /*!
      \brief estimates principal directions of the contour.

      This method estimates the principal contour directions automatically.

      \param input_range
      a range of input segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` 
      to `GeomTraits::Segment_2`

      \pre `input_range.size() > 1`
    */
    void estimate_principal_directions(
      const InputRange& input_range,
      const SegmentMap segment_map) {
      
      CGAL_precondition(input_range.size() > 1);
    }

    /*!
      \brief executes the contour regularization algorithm.

      This method regularizes the contour with respect to the defined principal directions.
      That is it sets all other not principal segments either orthogonal or collinear 
      to the chosen directions.

      \param input_range
      a range of input segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` 
      to `GeomTraits::Segment_2`

      \pre `input_range.size() > 1`
    */
    void regularize(
      InputRange& input_range,
      SegmentMap segment_map) {
      
      CGAL_precondition(input_range.size() > 1);
    }

  private:
    std::string m_parameters;
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
