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

// CGAL includes.
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// Eigen includes.
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
  */
  template<typename FT>
  class OSQP_quadratic_program : 
    public CGAL::Quadratic_program<int> {

    using Solution = CGAL::Quadratic_program_solution<FT>;
    using Triplet = Eigen::Triplet<FT>;
    using Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;
    using Sparse_matrix_iterator = typename Sparse_matrix::InnerIterator;

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

    void set_d(int i, int j, const FT& val) {
      P_vec.push_back(Triplet(i, j, val));
    }

    void set_c(int, const FT& val) {
      q_vec.push_back(val);
    }

    void set_c0(const FT&) {
      // empty!
    }

    void set_a(int j, int i, const FT& val) {
      A_vec.push_back(Triplet(i, j, val));
    }
    
    void set_b(int, const FT& val) {
      l_vec.push_back(-internal::max_value<FT>());
      u_vec.push_back(val);
    }

    void set_l(int, bool, const FT& val) {
      l_vec.push_back(val);
    }

    void set_u(int, bool, const FT& val) {
      u_vec.push_back(val);
    }

    bool solve(
      std::vector<FT>& solution) {

      finalize_qp_data();
      const c_int n = static_cast<c_int>(P_.nonZeros());
      const c_int m = static_cast<c_int>(l_.nonZeros());

      const c_int P_nnz = static_cast<c_int>(P_.nonZeros());
      const c_int q_nnz = static_cast<c_int>(q_.nonZeros());

      const c_int A_nnz = static_cast<c_int>(A_.nonZeros());
      const c_int l_nnz = static_cast<c_int>(l_.nonZeros());
      const c_int u_nnz = static_cast<c_int>(u_.nonZeros());

      CGAL_precondition(P_.nonZeros() == q_.nonZeros());
      CGAL_precondition(l_.nonZeros() == u_.nonZeros());

      c_float P_x[P_nnz];
      c_int   P_i[P_nnz];
      c_int   P_p[P_nnz + 1];
      set_P_data(P_nnz, P_, P_x, P_i, P_p);

      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[P_nnz + 1];
      set_A_data(A_, A_x, A_i, A_p);

      c_float q_x[q_nnz];
      c_float l_x[l_nnz];
      c_float u_x[u_nnz];
      set_qlu_data(n, m, q_, l_, u_, q_x, l_x, u_x);

      // Problem settings.
      OSQPSettings *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));

      // Structures.
      OSQPWorkspace *work;
      OSQPData *data;

      // Populate data.
      data = (OSQPData *)c_malloc(sizeof(OSQPData));
      data->n = n;
      data->m = m;
      data->P = csc_matrix(data->n, data->n, P_nnz, P_x, P_i, P_p);
      data->q = q_x;
      data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
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
      c_float *x = work->solution->x;
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
    std::vector<Triplet> P_vec;
    std::vector<Triplet> A_vec;
    std::vector<FT> q_vec, l_vec, u_vec;

    Sparse_matrix P_, A_;
    Dense_vector q_, l_, u_;

    void finalize_qp_data() {
      
      // P data.
      const std::size_t n = P_vec.size();
      P_.resize(n, n);
      P_.setFromTriplets(P_vec.begin(), P_vec.end());
      P_.makeCompressed();
      
      // A data.
      std::size_t s = A_vec.size() / 3;
      A_vec.reserve(A_vec.size() + n);
      for (std::size_t i = 0; i < n; ++i)
        A_vec.push_back(Triplet(s + i, i, FT(1)));

      const std::size_t m = s + n;
      A_.resize(m, n);
      A_.setFromTriplets(A_vec.begin(), A_vec.end());
      A_.makeCompressed(); 

      // q data.
      q_.resize(q_vec.size());
      for (std::size_t i = 0; i < q_vec.size(); ++i)
        q_[i] = q_vec[i];

      // l data.
      l_.resize(l_vec.size());
      for (std::size_t i = 0; i < l_vec.size(); ++i)
        l_[i] = l_vec[i];

      // u data.
      u_.resize(u_vec.size());
      for (std::size_t i = 0; i < u_vec.size(); ++i)
        u_[i] = u_vec[i];
    }

    void set_P_data(
      const c_int n, 
      const Sparse_matrix& P,
      c_float *P_x, 
      c_int   *P_i, 
      c_int   *P_p) const {

      std::size_t it = 0;
      for (int i = 0; i < P.outerSize(); ++i) {
        for (Sparse_matrix_iterator m_i(P, i); m_i; ++m_i) {
          const double val = CGAL::to_double(m_i.value());
          P_x[it] = val; ++it;
        }
      }
      P_p[0] = 0;
      for (int i = 0; i < n; ++i) {
        P_i[i] = i;
        P_p[i] = i;
      }
      P_p[n] = n;
    }

    void set_A_data(
      const Sparse_matrix& A,
      c_float *A_x, 
      c_int   *A_i, 
      c_int   *A_p) const {

      std::size_t it = 0;
      for (int i = 0; i < A.outerSize(); ++i) {
        for (Sparse_matrix_iterator m_i(A, i); m_i; ++m_i) {
          const double val = CGAL::to_double(m_i.value());
          const std::size_t idx = m_i.row();
          A_x[it] = val;
          A_i[it] = idx;
          ++it;
        }
      }
      A_p[0] = 0;
      for (int i = 1; i <= A.outerSize(); ++i) {
        const std::size_t coln = A.innerVector(i - 1).nonZeros();
        A_p[i] = A_p[i-1] + coln;
      } 
    }

    void set_qlu_data(
      const c_int n, const c_int m, 
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      c_float *q_x, 
      c_float *l_x, 
      c_float *u_x) const {

      for (int i = 0; i < m; ++i) {
        if (i < n) q_x[i] = CGAL::to_double(q[i]);
        l_x[i] = CGAL::to_double(l[i]);
        u_x[i] = CGAL::to_double(u[i]);
      }
    }
  };

  /*!
    \ingroup PkgShapeRegularizationRefSolvers
    
    \brief solves an OSQP quadratic program.

    \tparam FT
    number type.

    \param qp
    a quadratic program to be solved

    \param solution
    a vector with the solution
  */
  template<typename FT>
  bool solve_quadratic_program(
    CGAL::Shape_regularization::OSQP_quadratic_program<FT>& qp,
    std::vector<FT>& solution) {
    
    return qp.solve(solution);
  }

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_OSQP_QUADRATIC_PROGRAM_H
