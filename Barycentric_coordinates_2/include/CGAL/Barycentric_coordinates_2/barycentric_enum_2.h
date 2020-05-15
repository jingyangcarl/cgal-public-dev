// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_ENUM_2_H
#define CGAL_BARYCENTRIC_ENUM_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

namespace CGAL {

/*!
  \ingroup PkgBarycentricCoordinates2Ref

  The namespace `Barycentric_coordinates` contains implementations of all
  generalized barycentric coordinates: 2D, 3D, related enumerations, etc.
*/
namespace Barycentric_coordinates {

/// \name Computation Policies
/// @{

/*!
  `Computation_policy` provides a way to choose an asymptotic time complexity
  of the algorithm and its precision.
*/
enum class Computation_policy {

  /*!
    Computation is very precise but has a quadratic time complexity with respect
    to the number of the polygon vertices. In addition, we check a position of
    the query point with respect to the polygon and use different computation
    strategies for different positions. This is the default policy.
  */
  PRECISE_COMPUTATION_WITH_EDGE_CASES = 0,

  /*!
    Computation is very precise but has a quadratic time complexity with respect
    to the number of the polygon vertices. No extra checks are carried out.
  */
  PRECISE_COMPUTATION = 1,

  /*!
    Computation has a linear time complexity with respect to the number of the
    polygon vertices, but may suffer imprecision near the polygon boundary. In
    addition, we check a position of the query point with respect to the polygon
    and use different computation strategies for different positions.
  */
  FAST_COMPUTATION_WITH_EDGE_CASES = 2,

  /*!
    Computation has a linear time complexity with respect to the number of the
    polygon vertices, but may suffer imprecision near the polygon boundary.
    No extra checks are carried out.
  */
  FAST_COMPUTATION = 3,

  /*!
    The default policy is `PRECISE_COMPUTATION_WITH_EDGE_CASES`.
  */
  DEFAULT = PRECISE_COMPUTATION_WITH_EDGE_CASES
};

/// @}

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_ENUM_2_H
