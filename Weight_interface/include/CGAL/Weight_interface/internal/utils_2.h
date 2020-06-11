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

#ifndef CGAL_GENERALIZED_WEIGHTS_UTILS_2_H
#define CGAL_GENERALIZED_WEIGHTS_UTILS_2_H

// #include <CGAL/license/Weight_interface.h>

// STL includes.
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <iterator>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

namespace CGAL {
namespace Generalized_weights {
namespace internal {

template<typename GeomTraits>
const typename GeomTraits::FT area_3(
  const typename GeomTraits::Point_3& p,
  const typename GeomTraits::Point_3& q,
  const typename GeomTraits::Point_3& r,
  const GeomTraits traits) {

  return 0;
}

} // namespace internal
} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_WEIGHTS_UTILS_2_H
