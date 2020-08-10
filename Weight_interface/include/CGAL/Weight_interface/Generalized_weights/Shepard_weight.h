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

#ifndef CGAL_GENERALIZED_SHEPARD_WEIGHT_H
#define CGAL_GENERALIZED_SHEPARD_WEIGHT_H

// #include <CGAL/license/Weight_interface.h>

// Internal includes.
#include <CGAL/Weight_interface/internal/utils.h>

namespace CGAL {
namespace Generalized_weights {

  /*!
    \ingroup PkgWeightInterfaceRefWeights

    \brief Shepard weight.

    The full weight is computed as

    \f$w = \frac{1}{r^a}\f$

    with notations shown in the figure below and \f$a\f$ any real number
    being the power parameter.

    For \f$a = 1\f$ this weight is equal to the
    `CGAL::Generalized_weights::Inverse_distance_weight`.

    \cgalFigureBegin{shepard_weight, shepard.svg}
      Notation used for the Shepard weight.
    \cgalFigureEnd

    \tparam GeomTraits
    must be a model of `AnalyticTraits`.

    \cgalModels `AnalyticWeight_2`
  */
  template<typename GeomTraits>
  class Shepard_weight {

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
    Shepard_weight(
      const FT a = FT(1), // default is for inverse distance weight
      const GeomTraits traits = GeomTraits()) :
    m_p(a), m_traits(traits)
    { }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes the Shepard weight.
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
      \brief computes the Shepard weight.
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
    const FT m_p;
    const GeomTraits m_traits;

    const FT weight(
      const FT d) const {

      FT w = FT(0);
      CGAL_assertion(d != FT(0));
      if (d != FT(0)) {
        FT denom = d;
        if (m_p != FT(1))
          denom = internal::power(m_traits, d, m_p);
        w = FT(1) / denom;
      }
      return w;
    }
  };

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Shepard weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.
  */
  template<typename Point_2>
  decltype(auto) shepard_weight_2(
    const Point_2& q, const Point_2& t, const Point_2& r, const Point_2& p) {

    using Traits = typename Kernel_traits<Point_2>::Kernel;
    const Shepard_weight<Traits> shepard;
    return shepard(q, t, r, p);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Shepard weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param q
    a query point

    \param t
    the first neighbor

    \param r
    the second neighbor

    \param p
    the third neighbor

    \return the computed weight.
  */
  template<typename Point_3>
  decltype(auto) shepard_weight_3(
    const Point_3& q, const Point_3& t, const Point_3& r, const Point_3& p) {

    using Traits = typename Kernel_traits<Point_3>::Kernel;
    const Shepard_weight<Traits> shepard;
    return shepard(q, t, r, p);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Shepard weight for 2D points.

    This function infers a traits class `GeomTraits` from the `Point_2` type.

    \tparam Point_2
    must be `CGAL::Point_2<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.
  */
  template<typename Point_2>
  decltype(auto) shepard_weight_2(
    const Point_2& p, const Point_2& q) {

    Point_2 stub;
    return shepard_weight_2(p, stub, q, stub);
  }

  /*!
    \ingroup PkgWeightInterfaceRefFreeFunctions

    \brief computes the Shepard weight for 3D points.

    This function infers a traits class `GeomTraits` from the `Point_3` type.

    \tparam Point_3
    must be `CGAL::Point_3<GeomTraits>`.

    \param p
    the first point

    \param q
    the second point

    \return the computed weight.
  */
  template<typename Point_3>
  decltype(auto) shepard_weight_3(
    const Point_3& p, const Point_3& q) {

    Point_3 stub;
    return shepard_weight_3(p, stub, q, stub);
  }

} // namespace Generalized_weights
} // namespace CGAL

#endif // CGAL_GENERALIZED_SHEPARD_WEIGHT_H
