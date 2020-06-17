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

#ifndef CGAL_GENERALIZED_SHEPARD_WEIGHT_2_H
#define CGAL_GENERALIZED_SHEPARD_WEIGHT_2_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils_2.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRef2DWeights

    \brief 2D Shepard weight.

    The full weight is computed as

    \f$w = \frac{1}{r^a}\f$

    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 1\f$ this weight is equal to the
    `CGAL::Generalized_weights::Inverse_distance_weight_2`.

    \cgalFigureBegin{shepard_weight, shepard.svg}
      Notation used for the Shepard weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits_2`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Shepard_weight_2 {

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

      \param a
      the power parameter

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.
    */
    Shepard_weight_2(
      const FT a = FT(1), // default is for inverse distance weight
      const GeomTraits traits = GeomTraits()) :
    m_p(a), m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes 2D Shepard weight.
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
      \brief computes 2D Shepard weight.
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
    const FT m_p;
    const GeomTraits m_traits;

    const FT weight(
      const FT rj) const {

      FT w = FT(0);
      CGAL_assertion(rj != FT(0));
      if (rj != FT(0)) {
        FT denom = rj;
        if (m_p != FT(1))
          denom = internal::power(m_traits, rj, m_p);
        w = FT(1) / denom;
      }
      return w;
    }
  };

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_SHEPARD_WEIGHT_2_H
