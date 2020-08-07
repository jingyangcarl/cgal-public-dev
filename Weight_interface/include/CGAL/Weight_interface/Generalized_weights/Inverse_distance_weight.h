// Copyright (c) 2020 GeometryFactory SARL (France).
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

#ifndef CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_H
#define CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief Inverse distance weight.

    The full weight is computed as

    \f$w = \frac{1}{r}\f$

    with notations shown in the figure below. This weight is a special case of
    the `CGAL::Generalized_weights::Shepard_weight`.

    \cgalFigureBegin{inverse_distance_weight, inverse_distance.svg}
      Notation used for the inverse distance weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Inverse_distance_weight {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using GT = GeomTraits;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// 2D point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// 3D point type.
    typedef typename GeomTraits::Point_3 Point_3;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.
    */
    Inverse_distance_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the inverse distance weight.
    */
    const FT operator()(
      const Point_2& p,
      const Point_2&,
      const Point_2& q,
      const Point_2&) const {

      const FT d =
        internal::distance_2(m_traits, p, q);
      return weight(d);
    }

    /*!
      \brief computes the inverse distance weight.
    */
    const FT operator()(
      const Point_3& p,
      const Point_3&,
      const Point_3& q,
      const Point_3&) const {

      const FT d =
        internal::distance_3(m_traits, p, q);
      return weight(d);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    const FT weight(
      const FT d) const {

      FT w = FT(0);
      CGAL_assertion(d != FT(0));
      if (d != FT(0))
        w = FT(1) / d;
      return w;
    }
  };

  template<typename Point_2>
  decltype(auto) inverse_distance_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    Inverse_distance_weight<Traits> inverse_distance;
    return inverse_distance(q, t, r, p);
  }

  template<typename Point_3>
  decltype(auto) inverse_distance_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    Inverse_distance_weight<Traits> inverse_distance;
    return inverse_distance(q, t, r, p);
  }

  template<typename Point_2>
  decltype(auto) inverse_distance_weight_2(
    const Point_2& p, const Point_2& q) {

    Point_2 stub;
    return inverse_distance_weight_2(p, stub, q, stub);
  }

  template<typename Point_3>
  decltype(auto) inverse_distance_weight_3(
    const Point_3& p, const Point_3& q) {

    Point_3 stub;
    return inverse_distance_weight_3(p, stub, q, stub);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_H
