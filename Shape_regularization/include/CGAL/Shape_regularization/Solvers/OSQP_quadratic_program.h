// Copyright (c) 2019 GeometryFactory Sarl (France).
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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H
#define CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {  

  /*!
    \ingroup PkgShapeRegularizationRefSolvers
    
    \brief wraps the external OSQP solver.

    This class wraps the external \ref thirdpartyOSQP "OSQP solver" 
    and sets all its parameters to default.

    \tparam FT
    number type.
  */
  template<typename FT>
  class OSQP_quadratic_program : 
    public CGAL::Quadratic_program<int> {

    using Solution = CGAL::Quadratic_program_solution<FT>;

  public:
    
    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
    */
    OSQP_quadratic_program() 
    { }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void set_a(int j, int i, const FT& val) {

    }
    
    void set_b(int i, const FT& val) {
      
    }

    void set_l(int j, bool is_finite, const FT& val = FT(0)) {
      
    }

    void set_u(int j, bool is_finite, const FT& val = FT(0)) {
      
    }

    void set_c(int j, const FT& val) {
      
    }

    void set_c0(const FT& val) {
      
    }

    void set_d(int i, int j, const FT& val) {
      
    }

    Solution solve() const {
      return Solution();
    }
    /// \endcond
  };

  /*!
    \ingroup PkgShapeRegularizationRefSolvers
    
    \brief solves an OSQP quadratic program.

    \tparam FT
    number type.

    \param qp
    a quadratic program to be solved
  */
  template<typename FT>
  CGAL::Quadratic_program_solution<FT> solve_quadratic_program(
    const CGAL::Shape_regularization::OSQP_quadratic_program<FT>& qp) {
    return qp.solve();
  }

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H
