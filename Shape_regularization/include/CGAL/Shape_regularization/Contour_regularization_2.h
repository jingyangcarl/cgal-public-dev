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
  */
  template<
  typename GeomTraits,
  typename InputRange>
  class Contour_regularization_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param input_range 
      an input range with contour points

      \pre `input_range.size() > 1`
    */
    Contour_regularization_2(
      InputRange& input_range) :
    m_input_range(input_range) { 
      CGAL_precondition(input_range.size() > 1);
    }

    /// @}

    /// \name Access 
    /// @{

    /*!
      \brief executes the contour regularization algorithm.
    */
    void regularize() {

    }

  private:
    Input_range& m_input_range;
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_CONTOUR_REGULARIZATION_2_H
