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

#ifndef CGAL_GENERALIZED_VORONOI_REGION_WEIGHT_H
#define CGAL_GENERALIZED_VORONOI_REGION_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefRegions

    \brief Voronoi region weight.

    This weight is the area of the shaded region in the figure below. The region
    is formed by two midpoints of the edges incident to `q` and the circumcenter of
    the triangle `[p, q, r]`.

    \cgalFigureBegin{voronoi_region_weight, voronoi_cell.svg}
      Notation used for the Voronoi region weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.
  */
  template<typename GeomTraits>
  class Voronoi_region_weight {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using GT = GeomTraits;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

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
    Voronoi_region_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D Voronoi area weight.
    */
    template<typename Point_2>
    const FT operator()(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return weight_2(p, q, r);
    }

    /*!
      \brief computes 2D Voronoi area weight.
    */
    const FT operator()(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      return weight_3(p, q, r);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    template<typename Point_2>
    const FT weight_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      const auto circumcenter_2 =
        m_traits.construct_circumcenter_2_object();
      const Point_2 center =
        circumcenter_2(p, q, r);
      const Point_2 m1 =
        internal::barycenter_2(m_traits, q, r);
      const Point_2 m2 =
        internal::barycenter_2(m_traits, q, p);

      const FT A1 = internal::positive_area_2(m_traits, q, m1, center);
      const FT A2 = internal::positive_area_2(m_traits, q, center, m2);
      return weight(A1, A2);
    }

    const FT weight_3(
      const Point_3& p,
      const Point_3& q,
      const Point_3& r) const {

      const auto circumcenter_3 =
        m_traits.construct_circumcenter_3_object();
      const Point_3 center =
        circumcenter_3(p, q, r);
      const Point_3 m1 =
        internal::barycenter_3(m_traits, q, r);
      const Point_3 m2 =
        internal::barycenter_3(m_traits, q, p);

      const FT A1 = internal::positive_area_3(m_traits, q, m1, center);
      const FT A2 = internal::positive_area_3(m_traits, q, center, m2);
      return weight(A1, A2);
    }

    const FT weight(
      const FT A1, const FT A2) const {

      return A1 + A2;
    }
  };

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 2D triangle [p, q, r].

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.
  */
  template<typename Point_2>
  decltype(auto) voronoi_area_2(
    const Point_2& p, const Point_2& q, const Point_2& r) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    const Voronoi_region_weight<Traits> voronoi_area;
    return voronoi_area(p, q, r);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Voronoi area on a 3D triangle [p, q, r].

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \param r
    the third point

    \return the computed area.
  */
  template<typename Point_3>
  decltype(auto) voronoi_area_3(
    const Point_3& p, const Point_3& q, const Point_3& r) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    const Voronoi_region_weight<Traits> voronoi_area;
    return voronoi_area(p, q, r);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_VORONOI_REGION_WEIGHT_H
