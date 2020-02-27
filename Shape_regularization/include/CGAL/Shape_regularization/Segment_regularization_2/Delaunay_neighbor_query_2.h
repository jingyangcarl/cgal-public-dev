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

  /*!
    \ingroup PkgShapeRegularizationRef_2D

    \brief A neighbor query based on a Delaunay triangulation, which enables to 
    find the nearest neighbors in a set of `GeomTraits::Segment_2`.

    This class returns indices of the nearest neighbors of a query segment 
    in a set of 2D segments.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the `InputRange` 
    and value type is `GeomTraits::Segment_2`.

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

    using Delaunay_triangulation = DT;
    using Indices = std::vector<std::size_t>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` 
      to `GeomTraits::Segment_2`

      \pre `input_range.size() > 1`
    */
    Delaunay_neighbor_query_2(
      InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) { 
      
      CGAL_precondition(
        input_range.size() > 1);
      m_groups.resize(m_input_range.size());
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief implements `NeighborQuery::operator()()`.

      This operator returns indices of segments, which are direct neighbors of 
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

      \tparam ItemRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

      \tparam IndexMap 
      must be an `LvaluePropertyMap` whose key type is the value type of `ItemRange`
      and value type is `std::size_t`.

      \param item_range
      an instance of ItemRange

      \param index_map
      an instance of IndexMap that returns an index stored in the `item_range` 
      of the segment in the `input_range`

      \pre `item_range.size() > 1`
    */
    template<
    typename ItemRange, 
    typename IndexMap = CGAL::Identity_property_map<std::size_t> >
  	void add_group(
      const ItemRange& item_range, 
      const IndexMap index_map = IndexMap()) { 
      
      CGAL_precondition(item_range.size() > 1);
      if (item_range.size() < 2) return;
      
      Indices group;
      group.reserve(item_range.size());
      for (const auto& item : item_range) {
        const std::size_t seg_index = get(index_map, item);
        group.push_back(seg_index);
      }

      build_delaunay_triangulation(group);
      add_neighbors();
    }

    /// @}

    /// \name Internal data management
    /// @{ 

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_groups.clear();
      m_groups.resize(m_input_range.size());
    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;

    Delaunay_triangulation m_delaunay;
    std::vector<Indices> m_groups;

    void build_delaunay_triangulation(
      const Indices& group) {

      m_delaunay.clear();
      for (const std::size_t seg_index : group) {

        const auto& segment = get(
          m_segment_map, *(m_input_range.begin() + seg_index));
        const auto& source = segment.source();
        const auto& target = segment.target();
        
        const auto middle = internal::middle_point_2(source, target);
        auto vh = m_delaunay.insert(middle);
        vh->info() = seg_index;
      }
    }

    void add_neighbors() {

      CGAL_precondition(m_groups.size() == m_input_range.size());
      for (auto vit = m_delaunay.finite_vertices_begin(); 
      vit != m_delaunay.finite_vertices_end(); ++vit) {

        const std::size_t seg_index_1 = vit->info();
        CGAL_precondition(
          seg_index_1 >= 0 && seg_index_1 < m_groups.size());
        auto& neighbors = m_groups[seg_index_1];

        auto vc = m_delaunay.incident_vertices(vit);
        const auto end = vc;
        do {
          if (!m_delaunay.is_infinite(vc)) {
            const std::size_t seg_index_2 = vc->info();
            CGAL_precondition( 
              seg_index_2 >= 0 && seg_index_2 < m_input_range.size());
            neighbors.push_back(seg_index_2);
          }
          ++vc;
        } while (vc != end);
      }
      m_delaunay.clear();
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNAY_NEIGHBOR_QUERY_2_H
