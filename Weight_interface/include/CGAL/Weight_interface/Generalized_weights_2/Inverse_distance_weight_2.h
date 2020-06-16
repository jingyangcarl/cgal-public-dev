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

#ifndef CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_2_H
#define CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief 2D inverse distance weight.

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Inverse_distance_weight_2 {

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
    Inverse_distance_weight_2(
      const GeomTraits traits = GeomTraits()) :
    m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D inverse distance weight.
    */
    const FT operator()(
      const Point_2& query,
      const Point_2&,
      const Point_2& vj,
      const Point_2&) const {

      const FT rj =
        internal::distance_2(m_traits, query, vj);
      return weight(rj);
    }

    /*!
      \brief computes 2D inverse distance weight.
    */
    const FT operator()(
      const Point_3& query,
      const Point_3&,
      const Point_3& vj,
      const Point_3&) const {

      const FT rj =
        internal::distance_3(m_traits, query, vj);
      return weight(rj);
    }

    /// @}

  private:
    const GeomTraits m_traits;

    const FT weight(
      const FT rj) const {

      FT w = FT(0);
      CGAL_assertion(rj != FT(0));
      if (rj != FT(0))
        w = FT(1) / rj;
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_INVERSE_DISTANCE_WEIGHT_2_H
