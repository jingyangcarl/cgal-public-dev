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

#ifndef CGAL_GENERALIZED_TRIANGLE_CELL_WEIGHT_2_H
#define CGAL_GENERALIZED_TRIANGLE_CELL_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DAverage

    \brief 2D triangle cell weight.

    This weight is the area of the shaded triangle in the figure below.

    \cgalFigureBegin{triangle_cell_weight, triangle_cell.svg}
      Notation used for the triangle cell weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `HalfWeight_2`
  */
  template<typename GeomTraits>
  class Triangle_cell_weight_2 {

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
    Triangle_cell_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D triangle cell weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vj,
      const Point_2& vp) const {

      return weight_2(query, vj, vp);
    }

    /*!
      \brief computes 2D triangle cell weight.
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

      const auto area_2 =
        m_traits.compute_area_2_object();
      const FT Aj = area_2(vj, vp, query);
      return weight(Aj);
    }

    const FT weight_3(
      const Point_3& query,
      const Point_3& vj,
      const Point_3& vp) const {

      const FT Aj =
        internal::area_3(m_traits, vj, vp, query);
      return weight(Aj);
    }

    const FT weight(
      const FT Aj) const {

      const FT w = Aj;
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_TRIANGLE_CELL_WEIGHT_2_H
