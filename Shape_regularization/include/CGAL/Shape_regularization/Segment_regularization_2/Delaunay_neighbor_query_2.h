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

#ifndef CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2_H
#define CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2_H

// #include <CGAL/license/Shape_regularization.h>

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/property_map.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/utils.h>

namespace CGAL {
namespace Shape_regularization {

  /*!
    \ingroup PkgShapeRegularization2DReg

    \brief A neighbor query based on a Delaunay triangulation, which enables to 
    find the nearest neighbors in a set of `GeomTraits::Segment_2`.

    This class returns indices of the nearest neighbors of a query segment 
    in a set of 2D segments.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be an `LvaluePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Segment_2`.

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
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;

    using VB = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Traits>;
    using DS = CGAL::Triangulation_data_structure_2<VB>;
    using DT = CGAL::Delaunay_triangulation_2<Traits, DS>;

    using Delaunay_triangulation = DT;
    using Indices = std::vector<std::size_t>;
    using Vertex_circulator = typename DT::Vertex_circulator;
    using Indices_map = std::vector<Indices>;
    /// \endcond

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param segment_map
      an instance of `SegmentMap` that maps an item from `input_range` 
      to `Kernel::Segment_2`
    */
    Delaunay_neighbor_query_2(
      InputRange& input_range, 
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map) { 
      m_neighbor_map.resize(
        m_input_range.size());
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
      std::vector<std::size_t> & neighbors) { 

      neighbors.clear();
      if (m_input_range.size() <= 1) return;
      CGAL_precondition(
        m_neighbor_map.size() == m_input_range.size());      
      CGAL_precondition(
        query_index >= 0 && query_index < m_input_range.size());
      neighbors = m_neighbor_map[query_index];
    }

    /// @}

    /// \name Utilities
    /// @{ 

    /*!
      \brief adds a group of items with segments and then updates internal Delaunay 
      triangulation and neighbors of each inserted segment.

      \tparam ItemRange 
      must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

      \tparam IndexMap 
      must be an `LvaluePropertyMap` whose key type is the value type of item range
      and value type is `std::size_t`.

      \param item_range
      an instance of ItemRange

      \param index_map
      an instance of IndexMap that returns the index of a segment stored in the
      `item_range` with respect to the total `input_range`
    */
    template<
    typename ItemRange, 
    typename IndexMap = CGAL::Identity_property_map<std::size_t> >
  	void add_group(
      const ItemRange& item_range, 
      const IndexMap index_map = IndexMap()) { 
        
      Indices group;
      for (const auto& item : item_range) {
        const std::size_t seg_index = get(index_map, item);
        group.push_back(seg_index);
      }
      if (group.size() == 0) return;

      update_delaunay_triangulation(group);
      if (m_delaunay.number_of_vertices() > 1)
        update_neighbors();
    }

    /// @}

    /// \name Internal data management
    /// @{ 

    /*!
      \brief deletes the information about the neighbors for all the segments
    */

    void clear() {
      m_delaunay.clear();
      m_neighbor_map.clear();
      m_neighbor_map.resize(m_input_range.size());
    }

    /// @}

  private:
    const Input_range& m_input_range;
    const Segment_map m_segment_map;
    Delaunay_triangulation m_delaunay;
    Indices_map m_neighbor_map;

    void update_delaunay_triangulation(
      const Indices& group) {
      
      for (const std::size_t idx : group) {
        const auto& segment = get(
          m_segment_map, *(m_input_range.begin() + idx));
        const auto& source = segment.source();
        const auto& target = segment.target();
        const auto  middle = internal::middle_point_2(source, target);
        auto vh = m_delaunay.insert(middle);
        vh->info() = idx;
      }
    }

    void update_neighbors() {
      
      for (auto vit = m_delaunay.finite_vertices_begin(); 
      vit != m_delaunay.finite_vertices_end(); ++vit) {

        const std::size_t seg_index_1 = vit->info();
        CGAL_precondition(
          seg_index_1 >= 0 && seg_index_1 < m_input_range.size());
        auto& seg_indices = m_neighbor_map[seg_index_1];
        seg_indices.clear();

        auto vc = m_delaunay.incident_vertices(vit);
        const auto start = vc;
        do {
          if (!m_delaunay.is_infinite(vc)) {
            const std::size_t seg_index_2 = vc->info();
            CGAL_precondition( 
              seg_index_2 >= 0 && seg_index_2 < m_input_range.size());
            seg_indices.push_back(seg_index_2);
          }
          ++vc;
        } while (vc != start);
      } 
    }
  };

} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_DELAUNEY_NEIGHBOR_QUERY_2_H
