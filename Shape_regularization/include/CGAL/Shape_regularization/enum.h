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

#ifndef CGAL_SHAPE_REGULARIZATION_ENUM_H
#define CGAL_SHAPE_REGULARIZATION_ENUM_H

// #include <CGAL/license/Shape_regularization.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef
      
    \brief This label represents a type of estimated contour direction.
  */
  enum class Direction_type {
    
    /// The chosen direction is set to the longest segment of the input contour.
    LONGEST = 0,

    /// The chosen directions are estimated and may contain more than one direction.
    MULTIPLE = 1

  }; // Direction_type

} // Shape_regularization
} // CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ENUM_H
