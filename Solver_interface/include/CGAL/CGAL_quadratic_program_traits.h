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
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_QUADRATIC_PROGRAM_TRAITS_H
#define CGAL_QUADRATIC_PROGRAM_TRAITS_H

// CGAL includes.
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

namespace CGAL {

  /*!
    \ingroup PkgSolverInterfaceNLP

    \brief wraps the \cgal QP solver.

    This class wraps the internal \ref PkgQPSolver "CGAL solver"
    and sets all its parameters to default.

    \tparam FT
    number type.

    \cgalModels `QuadraticProgramTraits`
  */
  template<typename FT>
  class CGAL_quadratic_program_traits {
    using Quadratic_program = CGAL::Quadratic_program<FT>;

  public:

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
    */
    CGAL_quadratic_program_traits() :
      m_quadratic_program(
        CGAL::SMALLER,
        true,
        -internal::max_value<FT>(),
        true,
        +internal::max_value<FT>())
    { }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void reserve_d(const std::size_t) { }
    void reserve_c(const std::size_t) { }
    void reserve_a(const std::size_t) { }
    void reserve_b(const std::size_t) { }
    void reserve_l(const std::size_t) { }
    void reserve_u(const std::size_t) { }

    void set_d(
      const std::size_t i,
      const std::size_t j,
      const FT val) {
      m_quadratic_program.set_d(i, j, val);
    }

    void set_c(
      const std::size_t j,
      const FT val) {
      m_quadratic_program.set_c(j, val);
    }

    void set_c0(const FT val) {
      m_quadratic_program.set_c0(val);
    }

    void set_a(
      const std::size_t j,
      const std::size_t i,
      const FT val) {
      m_quadratic_program.set_a(j, i, val);
    }

    void set_b(
      const std::size_t i,
      const FT val) {
      m_quadratic_program.set_b(i, val);
    }

    void set_l(
      const std::size_t j,
      const bool is_finite,
      const FT val) {
      m_quadratic_program.set_l(j, is_finite, val);
    }

    void set_u(
      const std::size_t j,
      const bool is_finite,
      const FT val) {
      m_quadratic_program.set_u(j, is_finite, val);
    }

    /*
      \brief solves a quadratic program.

      \param solution
      a vector with the solution

      \returns a status of the computation `success == true`
    */
    bool solve(
      std::vector<FT>& solution) {

      auto s = CGAL::solve_quadratic_program(m_quadratic_program, FT());
      for (auto x = s.variable_values_begin();
      x != s.variable_values_end(); ++x)
        solution.push_back(static_cast<FT>(
          CGAL::to_double(*x)));

      return s.solves_quadratic_program(m_quadratic_program);
    }
    /// \endcond

  private:
    Quadratic_program m_quadratic_program;
  };

} // namespace CGAL

#endif // CGAL_QUADRATIC_PROGRAM_TRAITS_H
