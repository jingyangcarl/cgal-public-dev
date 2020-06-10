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

#include <CGAL/license/Barycentric_coordinates_2.h>
// #include <CGAL/license/Weight_interface.h>

namespace CGAL {

/*!
  \ingroup PkgWeightInterfaceRef

  The namespace `Generalized_weights` contains implementations of all
  generalized weights: 1D, 2D, 3D, related enumerations, etc.
*/
namespace Generalized_weights {

/// \name Computation Policies
/// @{

/*!
  `Computation_policy_2` provides a way to choose a type of algorithm and
  its precision for computing 2D generalized weights.
*/
enum class Computation_policy_2 {

  /*!
    The best trade-off linear-time algorithm for computing weights. In addition,
    we check a position of the query point with respect to the polygon and use
    different computation strategies for different positions. This is the default policy.
  */
  OPTIMAL_WITH_EDGE_CASES = 0,

  /*!
    The best trade-off linear-time algorithm for computing weights.
    No extra checks are carried out.
  */
  OPTIMAL = 1,

  /*!
    The default policy is `OPTIMAL_WITH_EDGE_CASES`.
  */
  DEFAULT = OPTIMAL_WITH_EDGE_CASES
};

/// @}

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_ENUM_2_H
