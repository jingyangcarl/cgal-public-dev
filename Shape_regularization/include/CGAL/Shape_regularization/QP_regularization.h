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
#include <Eigen/Sparse>
#include <Eigen/Dense>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/assertions.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularizationRef
    
    \brief Main class/entry point to the shape regularization algorithm
    based on the quadratic programming global optimization.

    Given a quadratic programming solver via `QPSolver`, this version of the 
    shape regularization algorithm enables to regularize a set of user-defined 
    items provided a way
    - to access neighbors of each item via the `NeighborQuery` class; 
    - to obtain a max bound for each item via the `RegularizationType` class;
    - to obtain a target value for each pair of neighbor items via the `RegularizationType` class;

    \tparam GeomTraits 
    must be a model of `Kernel`.
    
    \tparam InputRange 
    must be a model of `ConstRange`.

    \tparam NeighborQuery 
    must be a model of `NeighborQuery`.

    \tparam RegularizationType
    must be a model of `RegularizationType`.

    \tparam QPSolver
    must be a model of `QPSolver`.
  */
  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery, 
  typename RegularizationType,
  typename QPSolver>
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
    using QP_solver = QPSolver;

    using FT = typename Traits::FT;
    using Triplet = Eigen::Triplet<FT>;
    using Sparse_matrix = typename QP_solver::Sparse_matrix;
    using Dense_vector = typename QP_solver::Dense_vector;

    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.
      
      \param input_range 
      a range of input items for shape regularization

      \param neighbor_query 
      an instance of `NeighborQuery` that is used internally to 
      access item's neighbors

      \param regularization_type 
      an instance of `RegularizationType` that is used internally to 
      obtain bounds and target values of the items.

      \param qp_solver
      an instance of `QPSolver` to solve the quadratic programming problem.

      \pre `input_range.size() > 1`
    */
    QP_regularization(
      InputRange& input_range, 
      NeighborQuery& neighbor_query, 
      RegularizationType& regularization_type,
      QPSolver& qp_solver) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_regularization_type(regularization_type),
    m_qp_solver(qp_solver),
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

      m_graph.clear();
      build_graph_of_neighbors();
      CGAL_assertion(m_graph.size() > 0);
      if (m_graph.size() == 0) return;

      m_bounds.clear();
      m_bounds.reserve(m_input_range.size());
      m_max_bound = -FT(1);

      obtain_bounds();

      CGAL_assertion(m_max_bound > 0);
      CGAL_assertion(m_bounds.size() == m_input_range.size());

      m_targets.clear();

      obtain_targets();
      if (m_targets.size() == 0) return;
      CGAL_assertion(m_targets.size() > 0);

      build_OSQP_solver_data(); 

      std::vector<FT> result_qp;
      std::size_t n = m_input_range.size() + m_targets.size();
      result_qp.reserve(n);

      m_qp_solver.solve(
        m_input_range.size(), 
        m_targets.size(), 
        m_P_mat, m_A_mat, m_q, m_l, m_u, result_qp);
      CGAL_assertion(result_qp.size() == n);

      m_regularization_type.update(result_qp);
    }

    /// @}

  private:
    Input_range& m_input_range;
    Neighbor_query& m_neighbor_query;
    Regularization_type& m_regularization_type;
    QP_solver& m_qp_solver;

    std::set<Size_pair> m_graph;
    std::map<Size_pair, FT> m_targets;
    const Parameters m_parameters;
    
    // Variables for the OSQP solver:
    Sparse_matrix m_P_mat;
    Sparse_matrix m_A_mat;
    Dense_vector m_q;
    Dense_vector m_l;
    Dense_vector m_u;

    FT m_max_bound;
    std::vector<FT> m_bounds;

    void build_graph_of_neighbors() {
      
      Size_pair p;
      Indices neighbors;
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        neighbors.clear();
        m_neighbor_query(i, neighbors);
        for (const std::size_t index : neighbors) {
          i < index ? p = std::make_pair(i, index) : p = std::make_pair(index, i);
          m_graph.insert(p);
        }
      }
    }

    void obtain_bounds() {
      for (std::size_t i = 0; i < m_input_range.size(); ++i) {
        const FT bound = m_regularization_type.bound(i);
        CGAL_assertion(bound >= 0);
        m_bounds.push_back(bound);
        m_max_bound = CGAL::max(bound, m_max_bound);
      }
    }

    void obtain_targets() {
      for(const auto& pair : m_graph) {
        const std::size_t i = pair.first;
        const std::size_t j = pair.second; 

        const FT tar_val = m_regularization_type.target_value(i, j);
        if (CGAL::abs(tar_val) < 
          m_regularization_type.bound(i) + m_regularization_type.bound(j)) {
          m_targets[pair] = tar_val;
        }
      }
    }

    void build_quadratic_matrix(
      const std::size_t n, const std::size_t k) {
      
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

      m_P_mat.resize(n, n);
      m_P_mat.setFromTriplets(vec.begin(), vec.end());
      m_P_mat.makeCompressed();
    }

    void build_linear_part_vactor(
      const std::size_t n, const std::size_t k) {
      
      m_q.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        if (i < k) {
          m_q[i] = FT(0);
        } else {
          m_q[i] = m_parameters.lambda * m_parameters.weight / 
          (FT(4) * m_max_bound * FT(n - k));
        }
      }
    }

    void build_linear_constraints_matrix(
      const std::size_t n, 
      const std::size_t m, 
      const std::size_t k,
      const std::size_t e,
      const std::size_t A_nnz) {
      
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

      m_A_mat.resize(m, n);
      m_A_mat.setFromTriplets(vec.begin(), vec.end());
      m_A_mat.makeCompressed();
    }

    void build_bounds_vectors(
      const std::size_t m, 
      const std::size_t k, 
      const std::size_t e) {
      
      m_u.resize(m);
      m_l.resize(m);
      auto tit = m_targets.begin();

      for(std::size_t i = 0; i < m; ++i) {
        if (i < 2 * e) {
          const FT val = tit->second;
          if (i % 2 == 0) 
            m_u[i] = m_parameters.val_neg * val;
          else {
            m_u[i] = m_parameters.val_pos * val; ++tit;
          }
          m_l[i] = m_parameters.neg_inf;
        }
        else if (i < 2 * e + k) {
          m_l[i] = -1 * m_max_bound;
          m_u[i] = m_max_bound;
        } else {
          m_l[i] = m_parameters.neg_inf;
          m_u[i] = m_parameters.pos_inf;
        }
      }
    }

    void build_OSQP_solver_data() {
      const std::size_t k = m_input_range.size(); // k segments
      const std::size_t e = m_targets.size(); // e edges
      const std::size_t n = k + e; // number of variables
      const std::size_t m = 2 * e + n; // number of constraints
      const std::size_t A_nnz = 6 * e + n;  // number of entries in the constraint matrix

      build_quadratic_matrix(n, k);
      build_linear_part_vactor(n, k);
      build_linear_constraints_matrix(n, m, k, e, A_nnz);
      build_bounds_vectors(m, k, e);
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_QP_REGULARIZATION_H
