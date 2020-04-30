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

#ifndef CGAL_SHAPE_REGULARIZATION_DELAUNAY_NEIGHBOR_QUERY_2_H
#define CGAL_SHAPE_REGULARIZATION_DELAUNAY_NEIGHBOR_QUERY_2_H

// #include <CGAL/license/Shape_regularization.h>

// CGAL includes.
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief A neighbor query based on a Delaunay triangulation, which enables to 
    find the nearest neighbors in a set of `GeomTraits::Segment_2`.

    This class returns indices of the nearest neighbors of a query segment 
    in a set of 2D segments.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a model of `ReadablePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.

    \cgalModels `NeighborQuery`
  */
  template<
    typename GeomTraits, 
    typename InputRange, 
    typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Delaunay_neighbor_query_2 {

  public:
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Traits>;
    using DS = CGAL::Triangulation_data_structure_2<VB>;
    using DT = CGAL::Delaunay_triangulation_2<Traits, DS>;

    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    using Delaunay_triangulation = DT;
    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from input range to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 1`
    */
    Delaunay_neighbor_query_2(
      const InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) { 
      
      CGAL_precondition(input_range.size() > 1); clear();
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief returns indices of segments, which are direct neighbors of 
      the query segment.

      \param query_index
      an index of the query segment

      \param neighbors
      indices of segments, which are direct neighbors of the query segment

      \pre `query_index >= 0 && query_index < input_range.size()`
    */
    void operator()(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) { 
      
      neighbors.clear();
      CGAL_precondition(
        m_groups.size() == m_input_range.size());
      CGAL_precondition(
        query_index >= 0 && query_index < m_input_range.size());
      if (query_index >= m_input_range.size()) return;
      neighbors = m_groups[query_index];
    }

    /// @}

    /// \name Utilities
    /// @{ 

    /*!
      \brief inserts a group of segments from `input_range` and computes their 
      neighbors within the group.

      \tparam IndexRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.
      The value type is `std::size_t`.

      \param index_range
      a const range of segment indices

      \pre `index_range.size() > 1`
    */
    template<typename IndexRange>
  	void add_group(
      const IndexRange& index_range) { 
      if (index_range.size() < 2) 
        return;
      
      Indices group;
      group.reserve(index_range.size());
      for (const std::size_t seg_index : index_range)
        group.push_back(seg_index);

      build_delaunay_triangulation(group);
      add_neighbors();
      ++m_num_groups;
    }

    /*!
      \brief inserts all input segments from `input_range` as one unique group.

      For more details, 
      see `CGAL::Shape_regularization::Delaunay_neighbor_query_2::add_group()`.
    */
    void create_unique_group() {
      
      CGAL_precondition(m_input_range.size() > 1);
      if (m_input_range.size() < 2) return;

      Indices group(m_input_range.size());
      std::iota(group.begin(), group.end(), 0);
      add_group(group);
    }

    /// @}

    /// \name Internal data management
    /// @{ 

    /*!
      \brief clears all internal data structures.
    */
    void clear() {

      m_groups.clear();
      m_num_groups = 0;
      CGAL_precondition(m_input_range.size() > 1);
      if (m_input_range.size() < 2) return;
      m_groups.resize(m_input_range.size());
    }

    /// @}

    // EXTRA METHODS TO TEST THE CLASS!
    /// \cond SKIP_IN_MANUAL
    const std::size_t number_of_groups() const { 
      return m_num_groups;
    }

    const std::size_t number_of_neighbors() const {
      std::size_t num_neighbors = 0;
      for (const auto& group : m_groups)
        num_neighbors += group.size();
      return num_neighbors;
    }

    void get_edges(
      std::vector<Segment_2>& edges) {

      // Create a graph.
      std::set<Size_pair> graph;
      build_graph_of_neighbors(graph);

      // Set graph edges.
      edges.clear();
      edges.reserve(graph.size());
      for (const auto& edge : graph) {
        const auto p = get_middle_point(edge.first);
        const auto q = get_middle_point(edge.second);
        edges.push_back(Segment_2(p, q));
      }
    }
    /// \endcond

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;

    Delaunay_triangulation m_delaunay;
    std::vector<Indices> m_groups;
    std::size_t m_num_groups;

    void build_delaunay_triangulation(
      const Indices& group) {

      m_delaunay.clear();
      for (const std::size_t seg_index : group) {

        const auto& segment = get(
          m_segment_map, *(m_input_range.begin() + seg_index));
        const auto& source = segment.source();
        const auto& target = segment.target();
        const auto vh = m_delaunay.insert(
          internal::middle_point_2(source, target));
        vh->info() = seg_index;
      }
    }

    void add_neighbors() {

      CGAL_assertion(m_groups.size() == m_input_range.size());
      for (auto vit = m_delaunay.finite_vertices_begin(); 
      vit != m_delaunay.finite_vertices_end(); ++vit) {

        const std::size_t seg_index_1 = vit->info();
        CGAL_assertion(
          seg_index_1 >= 0 && seg_index_1 < m_groups.size());
        auto& neighbors = m_groups[seg_index_1];

        auto vc = m_delaunay.incident_vertices(vit);
        if (vc.is_empty()) return;
        const auto end = vc;
        do {
          if (!m_delaunay.is_infinite(vc)) {
            const std::size_t seg_index_2 = vc->info();
            CGAL_assertion( 
              seg_index_2 >= 0 && seg_index_2 < m_groups.size());
            neighbors.push_back(seg_index_2);
          }
          ++vc;
        } while (vc != end);
      }
      m_delaunay.clear();
    }

    // ALL BELOW USED ONLY FOR TESTING!
    void build_graph_of_neighbors(
      std::set<Size_pair>& graph) {
      
      Size_pair pair;
      Indices neighbors;
      graph.clear();

      CGAL_assertion(m_groups.size() == m_input_range.size());
      for (std::size_t i = 0; i < m_groups.size(); ++i) {
        operator()(i, neighbors);
        for (const std::size_t neighbor : neighbors) {
          i < neighbor ? 
          pair = std::make_pair(i, neighbor) : 
          pair = std::make_pair(neighbor, i);
          graph.insert(pair);
        }
      }
    }

    const Point_2 get_middle_point(
      const std::size_t seg_index) const {

      const auto& segment = get(
        m_segment_map, *(m_input_range.begin() + seg_index));
      const auto& source = segment.source();
      const auto& target = segment.target();
      return internal::middle_point_2(source, target);
    }
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNAY_NEIGHBOR_QUERY_2_H
