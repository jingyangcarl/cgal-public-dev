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

#ifndef CGAL_GENERALIZED_MIXED_VORONOI_REGION_WEIGHT_H
#define CGAL_GENERALIZED_MIXED_VORONOI_REGION_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DAverage

    \brief Mixed Voronoi region weight.

    This weight is the area of the shaded region in the figure below. The region
    is formed by two midpoints of the edges incident to `q` and the circumcenter of
    the triangle `[p, q, r]`.

    \cgalFigureBegin{mixed_voronoi_region_weight, mixed_voronoi_cell.svg}
      Notation used for the mixed Voronoi region weight.
    \cgalFigureEnd

    However, unlike the original `CGAL::Generalized_weights::Voronoi_region_weight`,
    if one of the angles in the triangle `[p, q, r]` is obtuse and the circumcenter
    vertex of the region is outside this triangle, this vertex is moved to the mid
    point of the edge `[r, p]`.

    \cgalFigureBegin{mixed_voronoi_region_weight_obtuse, mixed_voronoi_cell_obtuse.svg}
      The case with the obtuse angle.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.
  */
  template<typename GeomTraits>
  class Mixed_voronoi_region_weight {

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
    Mixed_voronoi_region_weight(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D mixed Voronoi area weight.
    */
    const FT operator()(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      return weight_2(p, q, r);
    }

    /*!
      \brief computes 2D mixed Voronoi area weight.
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

    const FT weight_2(
      const Point_2& p,
      const Point_2& q,
      const Point_2& r) const {

      const auto angle_2 =
        m_traits.angle_2_object();
      const auto a1 = angle_2(p, q, r);
      const auto a2 = angle_2(q, r, p);
      const auto a3 = angle_2(r, p, q);

      Point_2 center;
      if (a1 != CGAL::OBTUSE && a2 != CGAL::OBTUSE && a3 != CGAL::OBTUSE) {
        const auto circumcenter_2 =
          m_traits.construct_circumcenter_2_object();
        center = circumcenter_2(p, q, r);
      } else {
        center = internal::barycenter_2(m_traits, r, p);
      }

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

      const auto angle_3 =
        m_traits.angle_3_object();
      const auto a1 = angle_3(p, q, r);
      const auto a2 = angle_3(q, r, p);
      const auto a3 = angle_3(r, p, q);

      Point_3 center;
      if (a1 != CGAL::OBTUSE && a2 != CGAL::OBTUSE && a3 != CGAL::OBTUSE) {
        const auto circumcenter_3 =
          m_traits.construct_circumcenter_3_object();
        center = circumcenter_3(p, q, r);
      } else {
        center = internal::barycenter_3(m_traits, r, p);
      }

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

  template<typename Point_2>
  decltype(auto) mixed_voronoi_area_2(
    const Point_2& p, const Point_2& q, const Point_2& r) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    Mixed_voronoi_region_weight<Traits> mixed_voronoi_area;
    return mixed_voronoi_area(p, q, r);
  }

  template<typename Point_3>
  decltype(auto) mixed_voronoi_area_3(
    const Point_3& p, const Point_3& q, const Point_3& r) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    Mixed_voronoi_region_weight<Traits> mixed_voronoi_area;
    return mixed_voronoi_area(p, q, r);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_MIXED_VORONOI_REGION_WEIGHT_H
