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

#ifndef CGAL_SHAPE_REGULARIZATION_QP_REGULARIZATION_H
#define CGAL_SHAPE_REGULARIZATION_QP_REGULARIZATION_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <set>
#include <map>
#include <vector>
#include <utility>

// Eigen includes.
#include <Eigen/Dense>
#include <Eigen/Sparse>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>
#include <CGAL/Shape_regularization/Solvers/CGAL_quadratic_program.h>
#include <CGAL/Shape_regularization/Solvers/OSQP_quadratic_program.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef
    
    \brief Shape regularization algorithm based on the quadratic programming 
    global optimization.

    Given a quadratic programming solver via the class `QuadraticProgram`, this version of the 
    shape regularization algorithm enables to regularize a set of user-defined 
    items provided a way
    - to access neighbors of each item via the `NeighborQuery` class; 
    - to obtain a max bound for each item via the `RegularizationType` class;
    - to obtain a target value for each pair of neighbor items via the `RegularizationType` class.

    \tparam GeomTraits 
    must be a model of `Kernel`.
    
    \tparam InputRange 
    must be a model of `ConstRange`.

    \tparam NeighborQuery 
    must be a model of `NeighborQuery`.

    \tparam RegularizationType
    must be a model of `RegularizationType`.

    \tparam QuadraticProgram
    must be a model of `QuadraticProgram`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery, 
  typename RegularizationType,
  typename QuadraticProgram>
  class QP_regularization {

  private:
    class Parameters {
    public:
      using FT = typename GeomTraits::FT;

      const FT weight,  lambda;
      const FT neg_inf, pos_inf;
      const FT val_neg, val_pos;

      Parameters():
      weight(FT(100000)), 
      lambda(FT(4) / FT(5)), 
      neg_inf(-internal::max_value<FT>()),
      pos_inf(+internal::max_value<FT>()),
      val_pos(FT(+2) * lambda),
      val_neg(FT(-2) * lambda) 
      { }
    };

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Regularization_type = RegularizationType;
    using Quadratic_program = QuadraticProgram;

    using FT = typename Traits::FT;
    using Triplet = Eigen::Triplet<FT>;
    using Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>;
    using Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>;

    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param input_range 
      a const range of input items for shape regularization

      \param neighbor_query 
      an instance of `NeighborQuery` that is used internally to 
      access item's neighbors

      \param regularization_type 
      an instance of `RegularizationType` that is used internally to 
      obtain bounds and target values of the items

      \param quadratic_program
      an instance of `QuadraticProgram` to solve the quadratic programming problem

      \pre `input_range.size() > 1`
    */
    QP_regularization(
      const InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegularizationType& regularization_type,
      QuadraticProgram& quadratic_program) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_regularization_type(regularization_type),
    m_quadratic_program(quadratic_program),
    m_parameters(Parameters()),
    m_max_bound(-FT(1)) { 
      
      CGAL_precondition(input_range.size() > 1);
    }

    /// @}

    /// \name Access 
    /// @{

    /*!
      \brief runs the shape regularization algorithm.
    */
    void regularize() { 
      if (m_input_range.size() < 2) return;

      // Graph.
      m_graph.clear();
      build_graph_of_neighbors(m_graph);
      CGAL_assertion(m_graph.size() > 0);
      if (m_graph.size() == 0) return;

      // Bounds.
      m_bounds.clear();
      m_bounds.reserve(m_input_range.size());
      m_max_bound = -FT(1);
      obtain_bounds(m_bounds, m_max_bound);

      CGAL_assertion(m_max_bound > 0);
      CGAL_assertion(m_bounds.size() == m_input_range.size());

      // Targets.
      m_targets.clear();
      obtain_targets(m_targets);
      if (m_targets.size() == 0) return;
      CGAL_assertion(m_targets.size() > 0);

      // Solver data.
      build_solver_data(m_P, m_A, m_q, m_l, m_u, m_c0);

      // Solve.
      std::vector<FT> result_qp;
      std::size_t n = m_input_range.size() + m_targets.size();
      result_qp.reserve(n);

      solve_quadratic_program(  
        m_P, m_A, m_q, m_l, m_u, m_c0, 
        m_quadratic_program, result_qp);
      CGAL_assertion(result_qp.size() == n);

      // Update.
      m_regularization_type.update(result_qp);
    }

    /// @}

  private:
    const Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    Quadratic_program& m_quadratic_program;

    std::set<Size_pair> m_graph;
    std::map<Size_pair, FT> m_targets;
    const Parameters m_parameters;

    std::vector<FT> m_bounds;
    FT m_max_bound;

    Sparse_matrix m_P;
    Sparse_matrix m_A;
    Dense_vector m_q;
    Dense_vector m_l;
    Dense_vector m_u;
    FT m_c0;

    void build_graph_of_neighbors(
      std::set<Size_pair>& graph) {
      
      Size_pair p;
      Indices neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        for (const std::size_t index : neighbors) {
          i < index ? p = std::make_pair(i, index) : p = std::make_pair(index, i);
          graph.insert(p);
        }
      }
    }

    void obtain_bounds(
      std::vector<FT>& bounds,
      FT& max_bound) {
      
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const FT bound = m_regularization_type.bound(i);
        CGAL_assertion(bound >= 0);
        bounds.push_back(bound);
        max_bound = CGAL::max(bound, max_bound);
      }
    }

    void obtain_targets(
      std::map<Size_pair, FT>& targets) {
      
      for(const auto& pair : m_graph) {
        const std::size_t i = pair.first;
        const std::size_t j = pair.second; 

        const FT target = m_regularization_type.target(i, j);
        if (CGAL::abs(target) < 
          m_regularization_type.bound(i) + m_regularization_type.bound(j)) {
          targets[pair] = target;
        }
      }
    }

    void build_solver_data(
      Sparse_matrix& P, 
      Sparse_matrix& A,
      Dense_vector& q,
      Dense_vector& l,
      Dense_vector& u,
      FT& c0) const {
      
      const std::size_t k = m_input_range.size(); // k segments
      const std::size_t e = m_targets.size(); // e edges
      const std::size_t n = k + e; // number of variables
      const std::size_t m = 2 * e + n; // number of constraints
      const std::size_t A_nnz = 6 * e + n;  // number of entries in the constraint matrix

      build_quadratic_term(n, k, P);
      build_linear_term(n, k, q);
      build_constant_term(c0);

      build_constraint_matrix(n, m, k, e, A_nnz, A);
      build_constraint_bounds(m, k, e, l, u);
    }

    void build_quadratic_term(
      const std::size_t n, 
      const std::size_t k,
      Sparse_matrix& P) const {
      
      std::vector<Triplet> vec;
      vec.reserve(k);
      for(std::size_t i = 0; i < n; ++i) {
        FT val = FT(0);
        if (i < k) {
          val = FT(2) * m_parameters.weight * (FT(1) - m_parameters.lambda) / 
          (m_bounds[i] * m_bounds[i] * FT(k));
        }
        vec.push_back(Triplet(i, i, val));
      }
      CGAL_assertion(vec.size() == n);

      P.resize(n, n);
      P.setFromTriplets(vec.begin(), vec.end());
      P.makeCompressed();
    }

    void build_linear_term(
      const std::size_t n, 
      const std::size_t k,
      Dense_vector& q) const {
      
      q.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        if (i < k) {
          q[i] = FT(0);
        } else {
          q[i] = m_parameters.lambda * m_parameters.weight / 
          (FT(4) * m_max_bound * FT(n - k));
        }
      }
    }

    void build_constant_term(
      FT& c0) const {
      c0 = FT(0);
    }

    void build_constraint_matrix(
      const std::size_t n, 
      const std::size_t m, 
      const std::size_t k,
      const std::size_t e,
      const std::size_t A_nnz,
      Sparse_matrix& A) const {
      
      std::vector<Triplet> vec;
      vec.reserve(A_nnz);

      std::size_t it = 0;
      std::size_t ij = k;

      for (const auto& target : m_targets) {
        const std::size_t i = target.first.first;
        const std::size_t j = target.first.second;

        vec.push_back(Triplet(it, i, m_parameters.val_neg));
        vec.push_back(Triplet(it, j, m_parameters.val_pos));
        vec.push_back(Triplet(it, ij, -1));
        ++it;

        vec.push_back(Triplet(it, i, m_parameters.val_pos));
        vec.push_back(Triplet(it, j, m_parameters.val_neg));
        vec.push_back(Triplet(it, ij, -1));
        ++it; ++ij;
      }

      for (std::size_t i = e * 2; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          if (j == i - e * 2) {
            vec.push_back(Triplet(i, j, 1));
          }
        }
      }
      CGAL_assertion(vec.size() == A_nnz);

      A.resize(m, n);
      A.setFromTriplets(vec.begin(), vec.end());
      A.makeCompressed();
    }

    void build_constraint_bounds(
      const std::size_t m, 
      const std::size_t k, 
      const std::size_t e,
      Dense_vector& l,
      Dense_vector& u) const {
      
      u.resize(m);
      l.resize(m);
      auto tit = m_targets.begin();

      for(std::size_t i = 0; i < m; ++i) {
        if (i < 2 * e) {
          const FT val = tit->second;
          if (i % 2 == 0) 
            u[i] = m_parameters.val_neg * val;
          else {
            u[i] = m_parameters.val_pos * val; ++tit;
          }
          l[i] = m_parameters.neg_inf;
        }
        else if (i < 2 * e + k) {
          l[i] = -1 * m_max_bound;
          u[i] = m_max_bound;
        } else {
          l[i] = m_parameters.neg_inf;
          u[i] = m_parameters.pos_inf;
        }
      }
    }

    void solve_quadratic_program(
      const Sparse_matrix& D, 
      const Sparse_matrix& A,
      const Dense_vector& c,
      const Dense_vector& l,
      const Dense_vector& u,
      const FT c0,
      Quadratic_program& qp,
      std::vector<FT>& result) const {

      const std::size_t n = static_cast<std::size_t>(D.nonZeros());
      const std::size_t m = static_cast<std::size_t>(l.nonZeros());
      const std::size_t k = m - n;

      CGAL_precondition(D.nonZeros() == c.nonZeros());
      CGAL_precondition(l.nonZeros() == u.nonZeros());

      set_quadratic_term(D, qp);
      set_linear_term(c, qp);
      set_constant_term(c0, qp);

      set_constraint_matrix(k, A, qp);
      set_constraint_bounds(k, u, qp);

      set_bounds(k, l, u, qp);

      auto solution = CGAL::Shape_regularization::
        solve_quadratic_program(qp);

      result.clear();
      result.reserve(n);
      
      std::size_t i = 0;
      for (auto x = solution.variable_values_begin(); 
      x != solution.variable_values_end(); ++x, ++i)
        if (i < n) result.push_back(static_cast<FT>(
          CGAL::to_double(*x)));
    }

    void set_quadratic_term(
      const Sparse_matrix& D,
      Quadratic_program& qp) const {
      
      for (std::size_t i = 0; i < D.rows(); ++i)
        qp.set_d(i, i, D.coeff(i, i));
    }

    void set_linear_term(
      const Dense_vector& c,
      Quadratic_program& qp) const {
      
      for (std::size_t i = 0; i < c.rows(); ++i)
        qp.set_c(i, c[i]);
    }

    void set_constant_term(
      const FT c0,
      Quadratic_program& qp) const {
      
      qp.set_c0(c0);
    }

    void set_constraint_matrix(
      const std::size_t k, 
      const Sparse_matrix& A,
      Quadratic_program& qp) const {

      for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i < k) {
          for (std::size_t j = 0; j < A.cols(); ++j)
            qp.set_a(i, j, A.coeff(i, j));
        }
      }
    }

    void set_constraint_bounds(
      const std::size_t k,
      const Dense_vector& u,
      Quadratic_program& qp) const {

      for (std::size_t i = 0; i < u.rows(); ++i) {
        if (i < k)
          qp.set_b(i, u[i]);
      }
    }

    void set_bounds(
      const std::size_t k,
      const Dense_vector& l,
      const Dense_vector& u,
      Quadratic_program& qp) const {

      const std::size_t num = m_input_range.size();
      CGAL_assertion(l.rows() == u.rows());
      for (std::size_t i = 0; i < l.rows(); ++i) {
        if (i >= k) {
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

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_QP_REGULARIZATION_H
