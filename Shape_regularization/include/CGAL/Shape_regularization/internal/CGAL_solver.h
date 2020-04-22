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
// Author(s)     : Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_INTERNAL_CGAL_SOLVER_H
#define CGAL_SHAPE_REGULARIZATION_INTERNAL_CGAL_SOLVER_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>
#include <utility>

// Eigen includes.
#include <Eigen/Sparse>
#include <Eigen/Dense> 

// CGAL includes.
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// TODO:
// * Remove this file!

namespace CGAL {
namespace Shape_regularization {  
namespace internal {

  template<typename GeomTraits>
  class CGAL_solver { 

  public:
    using Traits = GeomTraits;
    
    using FT = typename GeomTraits::FT;
    using Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;
    using Sparse_matrix_iterator = typename Sparse_matrix::InnerIterator;
    using QP_problem = CGAL::Quadratic_program<int>;

    void solve(
      const std::size_t number_of_items,
      const std::size_t number_of_edges, 
      const Sparse_matrix& P, 
      const Sparse_matrix& A,
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      std::vector<FT>& result) {

      const std::size_t n = static_cast<std::size_t>(P.nonZeros());
      const std::size_t m = static_cast<std::size_t>(l.nonZeros());
      const std::size_t k = m - n;

      CGAL_precondition(P.nonZeros() == q.nonZeros());
      CGAL_precondition(l.nonZeros() == u.nonZeros());

      const FT neg_inf = -internal::max_value<FT>();
      const FT pos_inf = +internal::max_value<FT>();

      QP_problem qp(
        CGAL::SMALLER, true, neg_inf, true, pos_inf);
      build_P_data(
        P, qp);
      build_A_data(
        k, A, qp);
      build_vectors(
        k, number_of_items, q, l, u, qp);

      auto solution = CGAL::solve_quadratic_program(qp, FT());
      CGAL_assertion(solution.solves_quadratic_program(qp));

      result.clear();
      result.reserve(n);
      
      std::size_t i = 0;
      for (auto x = solution.variable_values_begin(); 
      x != solution.variable_values_end(); ++x, ++i)
        if (i < n) result.push_back(static_cast<FT>(
          CGAL::to_double(*x)));
    }

  private:
    void build_P_data(
      const Sparse_matrix& P,
      QP_problem& qp) const {

      for (std::size_t i = 0; i < P.rows(); ++i)
        qp.set_d(i, i, P.coeff(i, i));
    }

    void build_A_data(
      const std::size_t k, 
      const Sparse_matrix& A,
      QP_problem& qp) const {

      for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i < k) {
          for (std::size_t j = 0; j < A.cols(); ++j)
            qp.set_a(i, j, A.coeff(i, j));
        }
      }
    }

    void build_vectors(
      const std::size_t k,
      const std::size_t num,
      const Dense_vector& q,
      const Dense_vector& l,
      const Dense_vector& u,
      QP_problem& qp) const {

      for (std::size_t i = 0; i < q.rows(); ++i)
        qp.set_c(i, q[i]);
      qp.set_c0(0.0);

      CGAL_assertion(l.rows() == u.rows());
      for (std::size_t i = 0; i < l.rows(); ++i) {
        if (i < k) {
          qp.set_b(i, u[i]);
        } else {

          const std::size_t idx = i - k;
          if (idx < num) {
            qp.set_l(idx, true, l[i]);
            qp.set_u(idx, true, u[i]);
          } else {
            qp.set_l(idx, false, l[i]);
            qp.set_u(idx, false, u[i]);
          }
        }
      }
    }
  };

} // internal
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_CGAL_SOLVER_H
