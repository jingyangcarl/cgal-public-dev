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

#ifndef CGAL_GENERALIZED_WEIGHTS_ENUM_2_H
#define CGAL_GENERALIZED_WEIGHTS_ENUM_2_H

// #include <CGAL/license/Weight_interface.h>

namespace CGAL {

/*!
  \ingroup PkgWeightInterfaceRef

  The namespace `Generalized_weights` contains implementations of all
  generalized weights: 2D, 3D, related enumerations, etc.
*/
namespace Generalized_weights {

/// \name Computation Policies
/// @{

/*!
  `Computation_policy` provides a way to choose an asymptotic time complexity
  of the algorithm and its precision.
*/
enum class Computation_policy {

  /*!
    Computation has a linear time complexity with respect to the number of the
    polygon vertices, but may suffer imprecisions near the polygon boundary. In
    addition, we check a position of the query point with respect to the polygon
    and use different computation strategies for different positions.
  */
  FAST_COMPUTATION_WITH_EDGE_CASES = 0,

  /*!
    Computation has a linear time complexity with respect to the number of the
    polygon vertices, but may suffer imprecisions near the polygon boundary.
    No extra checks are carried out.
  */
  FAST_COMPUTATION = 1,

  /*!
    The default policy is `FAST_COMPUTATION_WITH_EDGE_CASES`.
  */
  DEFAULT = FAST_COMPUTATION_WITH_EDGE_CASES
};

/// @}

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_ENUM_2_H
