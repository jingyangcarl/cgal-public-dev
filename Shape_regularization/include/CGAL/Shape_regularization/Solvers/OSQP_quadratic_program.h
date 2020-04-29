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

#ifndef CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H
#define CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H

// #include <CGAL/license/Shape_regularization.h>

// OSQP includes.
#include <osqp/osqp.h>

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

    \cgalModels `QPSolver`
  */
  template<typename FT>
  class OSQP_quadratic_program {
    // row, col, value
    using Triplet = std::tuple<std::size_t, std::size_t, FT>;
    
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
    void reserve_d(const std::size_t n) {
      P_vec.reserve(n);
    }

    void reserve_c(const std::size_t n) {
      q_vec.reserve(n);
    }

    void reserve_a(const std::size_t n) {
      A_vec.reserve(n);
    }

    void reserve_b(const std::size_t n) {
      l_vec.reserve(n);
      u_vec.reserve(n);
    }

    void reserve_l(const std::size_t n) {
      l_vec.reserve(l_vec.size() + n);
    }

    void reserve_u(const std::size_t n) {
      u_vec.reserve(u_vec.size() + n);
    }

    void set_d(
      const std::size_t, 
      const std::size_t, 
      const FT val) {
      P_vec.push_back(val);
    }

    void set_c(
      const std::size_t, 
      const FT val) {
      q_vec.push_back(val);
    }

    void set_c0(const FT) {
      // It is not used by OSQP.
    }

    void set_a(
      const std::size_t j, 
      const std::size_t i, 
      const FT val) {
      A_vec.push_back(std::make_tuple(i, j, val));
    }
    
    void set_b(
      const std::size_t, 
      const FT val) {
      l_vec.push_back(-internal::max_value<FT>());
      u_vec.push_back(val);
    }

    void set_l(std::size_t, bool, const FT val) {
      l_vec.push_back(val);
    }

    void set_u(
      const std::size_t, 
      const bool, 
      const FT val) {
      u_vec.push_back(val);
    }

    /*
      \brief solves an OSQP quadratic program.

      \param solution
      a vector with the solution

      \returns a status of the computation `success == true`
    */
    bool solve(
      std::vector<FT>& solution) {

      finalize_qp_data();
      CGAL_precondition(P_vec.size() == q_vec.size());
      CGAL_precondition(l_vec.size() == u_vec.size());

      const c_int P_nnz = static_cast<c_int>(non_zeros_P());
      c_float P_x[P_nnz];
      c_int   P_i[P_nnz];
      c_int   P_p[P_nnz + 1];
      set_P_data(P_x, P_i, P_p);

      const c_int A_nnz = static_cast<c_int>(non_zeros_A());
      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[P_nnz + 1];
      set_A_data(A_x, A_i, A_p);

      const c_int q_nnz = static_cast<c_int>(non_zeros_q());
      const c_int l_nnz = static_cast<c_int>(non_zeros_l());
      const c_int u_nnz = static_cast<c_int>(non_zeros_u());

      c_float q_x[q_nnz];
      c_float l_x[l_nnz];
      c_float u_x[u_nnz];
      set_qlu_data(q_x, l_x, u_x);

      // Problem settings.
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures.
      OSQPWorkspace *work;
      OSQPData *data;

      // Populate data.
      const c_int n = P_nnz;
      const c_int m = l_nnz;

      data = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(n, n, P_nnz, P_x, P_i, P_p);
      data->q = q_x;
      data->A = csc_matrix(m, n, A_nnz, A_x, A_i, A_p);
      data->l = l_x;
      data->u = u_x;

      // Set solver settings.
      osqp_set_default_settings(settings);
      settings->eps_rel = 1.0e-15;
      settings->verbose = false;

      // Set workspace.
      osqp_setup(&work, data, settings);

      // Solve problem.
      const int exitflag = osqp_solve(work);
      const bool success = exitflag == 0 ? true : false;

      // Create solution.
      const c_float *x = work->solution->x;
      for (std::size_t i = 0; i < n; ++i) {
        const FT val = static_cast<FT>(x[i]);
        solution.push_back(val);
      }

      // Clean workspace.
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings);

      return success;
    }
    /// \endcond

  private:
    std::vector<Triplet> A_vec;
    std::vector<FT> P_vec, q_vec, l_vec, u_vec;

    void finalize_qp_data() {
      
      const std::size_t n = P_vec.size();
      const std::size_t s = A_vec.size() / 3;
      A_vec.reserve(A_vec.size() + n);
      for (std::size_t i = 0; i < n; ++i)
        A_vec.push_back(std::make_tuple(s + i, i, FT(1)));
    }

    const std::size_t non_zeros_P() const {
      return P_vec.size();
    }

    const std::size_t non_zeros_A() const {

      std::size_t count = 0;
      const std::size_t num_cols = std::get<1>(A_vec.back());
      for (std::size_t ref_col = 0; 
      ref_col <= num_cols; ++ref_col) {
        for (std::size_t i = 0; i < A_vec.size(); ++i) {
          const std::size_t col = std::get<1>(A_vec[i]);
          if (col == ref_col)
            ++count;
        }
      }
      return count;
    }

    const std::size_t non_zeros_q() const {
      return q_vec.size();
    }

    const std::size_t non_zeros_l() const {
      return l_vec.size();
    }

    const std::size_t non_zeros_u() const {
      return u_vec.size();
    }

    void set_P_data(
      c_float *P_x, 
      c_int   *P_i, 
      c_int   *P_p) const {

      const std::size_t n = P_vec.size();
      for (std::size_t i = 0; i < n; ++i)
        P_x[i] = CGAL::to_double(P_vec[i]);

      P_p[0] = 0;
      for (std::size_t i = 0; i < n; ++i) {
        P_i[i] = i;
        P_p[i] = i;
      }
      P_p[n] = n;
    }

    void set_A_data(
      c_float *A_x, 
      c_int   *A_i, 
      c_int   *A_p) const {

      A_p[0] = 0;
      std::size_t count = 0;
      const std::size_t num_cols = std::get<1>(A_vec.back());
      for (std::size_t ref_col = 0; 
      ref_col <= num_cols; ++ref_col) {
        
        std::size_t num_rows = 0;
        for (std::size_t i = 0; i < A_vec.size(); ++i) {
          const std::size_t row = std::get<0>(A_vec[i]);
          const std::size_t col = std::get<1>(A_vec[i]);
          const double val = CGAL::to_double(std::get<2>(A_vec[i]));

          if (col == ref_col) {
            A_i[count] = row; A_x[count] = val; 
            ++count; ++num_rows;
          }
        }
        A_p[ref_col + 1] = A_p[ref_col] + num_rows;
      }
    }

    void set_qlu_data(
      c_float *q_x, 
      c_float *l_x, 
      c_float *u_x) const {

      CGAL_assertion(l_vec.size() == u_vec.size());
      const std::size_t n = q_vec.size();
      const std::size_t m = l_vec.size();
      CGAL_assertion(n <= m);
      
      for (std::size_t i = 0; i < m; ++i) {
        if (i < n) q_x[i] = CGAL::to_double(q_vec[i]);
        l_x[i] = CGAL::to_double(l_vec[i]);
        u_x[i] = CGAL::to_double(u_vec[i]);
      }
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H
