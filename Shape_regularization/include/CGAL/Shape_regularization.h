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
// Author(s)     : Dmitry Anisimov, Gennadii Sytov, Simon Giraudot, Jean-Philippe Bauchet, and Florent Lafarge
//

#ifndef CGAL_SHAPE_REGULARIZATION_HEADERS_H
#define CGAL_SHAPE_REGULARIZATION_HEADERS_H

// #include <CGAL/license/Shape_regularization.h>

#include <CGAL/Shape_regularization/QP_regularization.h>
#include <CGAL/Shape_regularization/Contour_regularization_2.h>
#include <CGAL/Shape_regularization/regularize_planes.h>

#include <CGAL/Shape_regularization/Segments/Delaunay_neighbor_query_2.h> 
#include <CGAL/Shape_regularization/Segments/Angle_regularization_2.h>
#include <CGAL/Shape_regularization/Segments/Offset_regularization_2.h>
#include <CGAL/Shape_regularization/Segments/Parallel_groups_2.h>

#include <CGAL/Shape_regularization/Contours/Longest_direction_2.h>
#include <CGAL/Shape_regularization/Contours/Multiple_directions_2.h>
#include <CGAL/Shape_regularization/Contours/User_defined_directions_2.h>

#include <CGAL/Shape_regularization/Solvers/CGAL_quadratic_program.h>
#include <CGAL/Shape_regularization/Solvers/OSQP_quadratic_program.h>

#endif // CGAL_SHAPE_REGULARIZATION_HEADERS_H
