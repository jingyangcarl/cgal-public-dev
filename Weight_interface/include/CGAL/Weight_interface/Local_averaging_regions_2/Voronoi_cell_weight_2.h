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

#ifndef CGAL_GENERALIZED_VORONOI_CELL_WEIGHT_2_H
#define CGAL_GENERALIZED_VORONOI_CELL_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DAverage

    \brief 2D Voronoi cell weight.

    This weight is the area of the shaded region in the Figure below. The region
    is formed by two mid points of the edges incident to `q` and the circumcenter of
    the triangle `[vj, vp, q]`.

    \cgalFigureBegin{voronoi_cell_weight, voronoi_cell.svg}
      Notation used for the Voronoi cell weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `HalfWeight_2`
  */
  template<typename GeomTraits>
  class Voronoi_cell_weight_2 {

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
    Voronoi_cell_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D Voronoi cell weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vj, vp);
    }

    /*!
      \brief computes 2D Voronoi cell weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      return weight_3(query, vj, vp);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    const FT weight_2(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      const auto circumcenter_2 =
        m_traits.construct_circumcenter_2_object();
      const Point_2 center =
        circumcenter_2(vj, vp, query);
      const Point_2 m1 =
        internal::barycenter_2(m_traits, query, vj);
      const Point_2 m2 =
        internal::barycenter_2(m_traits, query, vp);

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT A1 = area_2(m1, center, query);
      const FT A2 = area_2(center, m2, query);
      return weight(A1, A2);
    }

    const FT weight_3(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      const auto circumcenter_3 =
        m_traits.construct_circumcenter_3_object();
      const Point_3 center =
        circumcenter_3(vj, vp, query);
      const Point_3 m1 =
        internal::barycenter_3(m_traits, query, vj);
      const Point_3 m2 =
        internal::barycenter_3(m_traits, query, vp);

      const FT A1 = internal::area_3(m_traits, m1, center, query);
      const FT A2 = internal::area_3(m_traits, center, m2, query);
      return weight(A1, A2);
    }

    const FT weight(
      const FT A1, const FT A2) const {

      const FT w = A1 + A2;
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_VORONOI_CELL_WEIGHT_2_H
