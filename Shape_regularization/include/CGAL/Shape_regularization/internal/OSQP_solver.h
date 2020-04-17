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
// Author(s)     : Gennadii Sytov, Dmitry Anisimov
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

// TODO:
// * Remove this file!

namespace CGAL {
namespace Shape_regularization {
namespace internal {

  template<typename GeomTraits>
  class OSQP_solver { 

  public:
    using Traits = GeomTraits;
    
    using FT = typename GeomTraits::FT;
    using Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;
    using Sparse_matrix_iterator = typename Sparse_matrix::InnerIterator;

    void solve(
      const std::size_t number_of_items,
      const std::size_t number_of_edges, 
      const Sparse_matrix& P, 
      const Sparse_matrix& A,
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      std::vector<FT>& result) {

      const c_int n = static_cast<c_int>(P.nonZeros());
      const c_int m = static_cast<c_int>(l.nonZeros());

      const c_int P_nnz = static_cast<c_int>(P.nonZeros());
      const c_int q_nnz = static_cast<c_int>(q.nonZeros());

      const c_int A_nnz = static_cast<c_int>(A.nonZeros());
      const c_int l_nnz = static_cast<c_int>(l.nonZeros());
      const c_int u_nnz = static_cast<c_int>(u.nonZeros());

      CGAL_precondition(P.nonZeros() == q.nonZeros());
      CGAL_precondition(l.nonZeros() == u.nonZeros());

      c_float P_x[P_nnz];
      c_int   P_i[P_nnz];
      c_int   P_p[P_nnz + 1];
      build_P_data(P_nnz, P, P_x, P_i, P_p);

      c_float A_x[A_nnz];
      c_int   A_i[A_nnz];
      c_int   A_p[P_nnz + 1];
      build_A_data(A, A_x, A_i, A_p);

      c_float q_x[q_nnz];
      c_float l_x[l_nnz];
      c_float u_x[u_nnz];
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
      result.reserve(n);
      
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

  private:
    void build_P_data(
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

    void build_A_data(
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
      c_float *u_x) const {

      for (int i = 0; i < m; ++i) {
        if (i < n) q_x[i] = CGAL::to_double(q[i]);
        l_x[i] = CGAL::to_double(l[i]);
        u_x[i] = CGAL::to_double(u[i]);
      }
    }
  };

} // namespace internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_OSQP_SOLVER_H
