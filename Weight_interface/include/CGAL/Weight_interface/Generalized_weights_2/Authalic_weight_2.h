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

#ifndef CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H
#define CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Boost includes.
#include <CGAL/boost/graph/helpers.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Wachspress_cot_weight_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2D

    \brief 2D authalic weight.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Authalic_weight_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using GT = GeomTraits;
    using Base = Wachspress_cot_weight_2<GeomTraits>;
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
    Authalic_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_base(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D authalic weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2& vm,
      const Point_2& vj,
      const Point_2& vp) const {

      return m_base(query, vm, vj, vp);
    }

    /*!
      \brief computes 2D authalic weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3& vm,
      const Point_3& vj,
      const Point_3& vp) const {

      return m_base(query, vm, vj, vp);
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    template<
    typename PolygonMesh,
    typename VertexDescriptor,
    typename VertextAroundTargetCirculator>
    const FT operator()(
      const PolygonMesh& polygon_mesh,
      const VertexDescriptor vdi,
      const VertextAroundTargetCirculator vcj) const {

      return m_base(polygon_mesh, vdi, vcj);
    }
    /// \endcond

  private:
    const Base m_base;
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_AUTHALIC_WEIGHT_2_H
