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

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// Eigen includes.
#include <Eigen/Sparse>
#include <Eigen/Dense> 

// OSQP includes.
#include <osqp/osqp.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularization
    
    \brief Quadratic programming solver

    This is the CGAL wrapper for the OSQP library: https://osqp.org/

    It is recommended to use this solver instead of the default
    `CGAL::Shape_regularization::CGAL_solver`.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \cgalModels `QPSolver`
  */
  template<typename GeomTraits>
  class OSQP_solver { 

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    /// \endcond
    
    // Number type.
    using FT = typename GeomTraits::FT;

    // Sparse matrix type.
    using Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;

    // Dense vector type.
    using Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;

    /// \cond SKIP_IN_MANUAL
    using Sparse_matrix_iterator = typename Sparse_matrix::InnerIterator;
    /// \endcond

    /// @}

    /// \name Solver
    /// @{ 

    /*!
      \brief implements `QPSolver::solve()`.

      This function computes the quadratic programming problem given its
      input data.

      \param number_of_items
      number of input items

      \param number_of_edges
      number of edges in the connectivity graph

      \param P
      quadratic term

      \param A
      linear constraints

      \param q
      linear term

      \param l
      lower bounds

      \param u
      upper bounds

      \pre `P.nonZeros() == number_of_items + number_of_edges`
      \pre `A.nonZeros() == 6 * number_of_edges + n`
      \pre `q.nonZeros() == number_of_items + number_of_edges`
      \pre `l.nonZeros() == 2 * number_of_edges + n`
      \pre `u.nonZeros() == 2 * number_of_edges + n`
    */
    void solve(
      const std::size_t number_of_items,
      const std::size_t number_of_edges, 
      const Sparse_matrix& P, 
      const Sparse_matrix& A,
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      std::vector<FT>& result) {

      const c_int n = number_of_items + number_of_edges; // number of variables
      const c_int m = 2 * number_of_edges + n; // number of constraints
      const c_int P_nnz = n;
      const c_int A_nnz = 6 * number_of_edges + n;

      CGAL_precondition(P.nonZeros() == n);
      CGAL_precondition(A.nonZeros() == A_nnz);
      CGAL_precondition(q.nonZeros() == n);
      CGAL_precondition(l.nonZeros() == m);
      CGAL_precondition(u.nonZeros() == m);

      c_float P_x[n];
      c_int   P_i[n];
      c_int   P_p[n+1];
      build_P_data(n, P, P_x, P_i, P_p);

      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[n+1];
      build_A_data(A, A_x, A_i, A_p);

      c_float q_x[n];
      c_float l_x[m];
      c_float u_x[m];
      build_vectors(n, m, q, l, u, q_x, l_x, u_x);

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

      // Define solver settings as default.
      osqp_set_default_settings(settings);
      settings->eps_rel = 1.0e-15;
      settings->verbose = false;

      // Setup workspace.
      const c_int exitflag = osqp_setup(&work, data, settings);

      // Solve problem.
      osqp_solve(work);

      result.clear();
      c_float *x = work->solution->x;
      for(int i = 0; i < n; ++i)
        result.push_back(static_cast<FT>(x[i]));

      // Clean workspace.
      osqp_cleanup(work);
      c_free(data->A);
      c_free(data->P);
      c_free(data);
      c_free(settings);
    }

    /// @}

  private:
    void build_P_data(
      const c_int n, 
      const Sparse_matrix& P,
      c_float *P_x, 
      c_int   *P_i, 
      c_int   *P_p) {

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

    void build_A_data(
      const Sparse_matrix& A,
      c_float *A_x, 
      c_int   *A_i, 
      c_int   *A_p) {

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
        const std::size_t coln = A.innerVector(i-1).nonZeros();
        A_p[i] = A_p[i-1] + coln;
      } 
    }

    void build_vectors(
      const c_int n, const c_int m, 
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      c_float *q_x, 
      c_float *l_x, 
      c_float *u_x) {

      for (int i = 0; i < m; ++i) {
        if (i < n) q_x[i] = CGAL::to_double(q[i]);
        l_x[i] = CGAL::to_double(l[i]);
        u_x[i] = CGAL::to_double(u[i]);
      }
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER_H
