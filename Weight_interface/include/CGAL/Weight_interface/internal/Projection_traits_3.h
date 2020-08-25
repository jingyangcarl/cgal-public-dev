// Copyright (c) 2020 GeometryFactory SARL (France)
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

#ifndef CGAL_WEIGHT_INTERFACE_INTERNAL_PROJECTION_TRAITS_3_H
#define CGAL_WEIGHT_INTERFACE_INTERNAL_PROJECTION_TRAITS_3_H

#include <CGAL/Filtered_predicate_with_state.h>
#include <CGAL/Weight_interface/internal/Projection_traits_base_3.h>

namespace CGAL {
namespace Weights {
namespace internal {

template<class Filtered_kernel >
class Filtered_projection_traits_3
  : public Projection_traits_base_3<Filtered_kernel> {

  using K    = Filtered_kernel;
  using Self = Filtered_projection_traits_3<K>;
  using Base = Projection_traits_base_3<K>;

  using Exact_kernel       = typename K::Exact_kernel;
  using Approximate_kernel = typename K::Approximate_kernel;

  using C2E = typename K::C2E;
  using C2F = typename K::C2F;

public:
  using Exact_traits     = Projection_traits_base_3<Exact_kernel>;
  using Filtering_traits = Projection_traits_base_3<Approximate_kernel>;

public:
  explicit Filtered_projection_traits_3(
    const typename K::Vector_3& normal) :
  Base(normal)
  { }

#define CGAL_PROJ_TRAITS_FILTER_PRED(P, Pf, ACCESSOR) \
  typedef Filtered_predicate_with_state< \
    typename Exact_traits::P, \
    typename Filtering_traits::P, \
    C2E, \
    C2F, \
    typename K::Vector_3> P; \
  P Pf() const { \
    return P(this->ACCESSOR()); \
  }

  CGAL_PROJ_TRAITS_FILTER_PRED(
    Orientation_2,
    orientation_2_object,
    normal)
  CGAL_PROJ_TRAITS_FILTER_PRED(
    Collinear_2,
    collinear_2_object,
    normal)
};

// This declaration is necessary for breaking the cyclic dependency.
template<class Filtered_kernel>
class Filtered_projection_traits_3;

template<class Kernel, bool Has_filtered_predicates = Kernel::Has_filtered_predicates>
class Projection_traits_3
  : public Projection_traits_base_3<Kernel> {
public:
  explicit
  Projection_traits_3(
    const typename Kernel::Vector_3& normal)
    : Projection_traits_base_3<Kernel>(normal)
  { }
};

template<class Kernel>
class Projection_traits_3<Kernel, true>
  : public Filtered_projection_traits_3<Kernel> {
public:
  explicit
  Projection_traits_3(
    const typename Kernel::Vector_3& normal) :
  Filtered_projection_traits_3<Kernel>(normal)
  { }
};

} // namespace internal
} // namespace Weights
} // namespace CGAL

#endif // CGAL_WEIGHT_INTERFACE_INTERNAL_PROJECTION_TRAITS_3_H
