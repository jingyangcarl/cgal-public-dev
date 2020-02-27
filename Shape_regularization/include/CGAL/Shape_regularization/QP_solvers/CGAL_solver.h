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
//

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef_Solvers
    
    \brief Quadratic programming solver

    This is the wrapper for the CGAL QP solver. This solver is the default,
    but slow. It is recommended to use the `CGAL::Shape_regularization::OSQP_solver`.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \cgalModels `QPSolver`
  */
  template<typename GeomTraits>
  class CGAL_solver { 

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

      \param result
      stores the optimization results

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

      const std::size_t n = number_of_items + number_of_edges; // number of variables
      const std::size_t m = 2 * number_of_edges + n; // number of constraints
      const std::size_t P_nnz = n;
      const std::size_t A_nnz = 6 * number_of_edges + n;

      CGAL_precondition(P.nonZeros() == n);
      CGAL_precondition(A.nonZeros() == A_nnz);
      CGAL_precondition(q.nonZeros() == n);
      CGAL_precondition(l.nonZeros() == m);
      CGAL_precondition(u.nonZeros() == m);

    }

  private:

  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_INTERNAL_CGAL_SOLVER_H
