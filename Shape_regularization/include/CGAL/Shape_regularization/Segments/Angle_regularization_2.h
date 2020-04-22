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
// Author(s)     : Jean-Philippe Bauchet, Florent Lafarge, Gennadii Sytov, Dmitry Anisimov
//

#ifndef CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H
#define CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H

// #include <CGAL/license/Shape_regularization.h>

// Boost includes.
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

// Internal includes.
#include <CGAL/Shape_regularization/internal/Segment_data_2.h>
#include <CGAL/Shape_regularization/internal/Grouping_segments_2.h>
#include <CGAL/Shape_regularization/internal/Angle_conditions_2.h>

// TODO:
// * Clean it up.

namespace CGAL {
namespace Shape_regularization {
namespace Segments {

  /*!
    \ingroup PkgShapeRegularizationRefSegments

    \brief An angle-based regularization type for 2D segments that preserves 
    parallelism and orthogonality relationships.

    \tparam GeomTraits 
    must be a model of `Kernel`.

    \tparam InputRange 
    must be a model of `ConstRange` whose iterator type is `RandomAccessIterator`.

    \tparam SegmentMap 
    must be a `ReadablePropertyMap` whose key type is the value type of the input 
    range and value type is `GeomTraits::Segment_2`. %Default is the 
    `CGAL::Identity_property_map<typename GeomTraits::Segment_2>`.

    \cgalModels `RegularizationType`
  */
  template<
  typename GeomTraits, 
  typename InputRange,
  typename SegmentMap = CGAL::Identity_property_map<typename GeomTraits::Segment_2> >
  class Angle_regularization_2 {
  public:

    /// \name Types
    /// @{
  
    /// \cond SKIP_IN_MANUAL
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Segment_map = SegmentMap;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// \cond SKIP_IN_MANUAL
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Segment_2 = typename Traits::Segment_2;

    using Segment_data = typename internal::Segment_data_2<Traits>;
    using Conditions = typename internal::Angle_conditions_2<Traits>;
    using Grouping = internal::Grouping_segments_2<Traits, Conditions>;
    
    using Indices = std::vector<std::size_t>;
    using Size_pair = std::pair<std::size_t, std::size_t>;
    
    using Targets_map = 
      std::map<Size_pair, std::pair< FT, std::size_t> >;
    using Relations_map = 
      std::map<Size_pair, std::pair<int, std::size_t> >;
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \tparam NamedParameters
      a sequence of \ref pmp_namedparameters "Named Parameters".

      \param input_range 
      an instance of `InputRange` with 2D segments

      \param np
      optional sequence of \ref pmp_namedparameters "Named Parameters" 
      among the ones listed below

      \param max_angle
      max angle bound in degrees, the default is 25 degrees 

      \param segment_map
      an instance of `SegmentMap` that maps an item from input range to `GeomTraits::Segment_2`, 
      if not provided, the default is used

      \pre `input_range.size() > 1`
      \pre `max_angle >= 0 && max_angle <= 90`
    */
    template<typename NamedParameters>
    Angle_regularization_2(
      InputRange& input_range, 
      const NamedParameters np,
      const SegmentMap segment_map = SegmentMap()) :
    m_input_range(input_range),
    m_segment_map(segment_map),
    m_num_modified_segments(0) { 
      
      FT max_angle = parameters::choose_parameter(
        parameters::get_parameter(np, internal_np::max_angle), FT(25));
      m_theta_max = max_angle;

      CGAL_precondition(input_range.size() > 1);
      CGAL_precondition(max_angle >= FT(0) && max_angle <= FT(90));

      if (max_angle < FT(0) || max_angle > FT(90)) {
        std::cout << 
          "WARNING: The max angle bound has to be within [0, 90]! Setting to 0." 
        << std::endl;
        m_theta_max = FT(0);
      }
    }

    /// @}

    /// \name Access
    /// @{ 

    /*!
      \brief calculates the target value between 2 segments, which are
      direct neighbors to each other. The traget value is the angle.

      \param query_index_i
      index of the first segment

      \param query_index_j
      index of the second segment

      \pre `query_index_i >= 0 && query_index_i < input_range.size()`
      \pre `query_index_j >= 0 && query_index_j < input_range.size()`
    */
    FT target(
      const std::size_t query_index_i, 
      const std::size_t query_index_j) {

      CGAL_precondition(m_segments.size() > 1);
      CGAL_precondition(m_segments.find(query_index_i) != m_segments.end());
      CGAL_precondition(m_segments.find(query_index_j) != m_segments.end());

      const std::size_t i = query_index_i;
      const std::size_t j = query_index_j;

      const auto& s_i = m_segments.at(i);
      const auto& s_j = m_segments.at(j);

      const FT mes_ij = s_i.orientation - s_j.orientation;
      const double mes_90 = std::floor(CGAL::to_double(mes_ij / FT(90)));

      const FT to_lower = FT(90) *  static_cast<FT>(mes_90)          - mes_ij;
      const FT to_upper = FT(90) * (static_cast<FT>(mes_90) + FT(1)) - mes_ij;

      const FT tar_val = 
        CGAL::abs(to_lower) < CGAL::abs(to_upper) ? to_lower : to_upper;

      if (CGAL::abs(tar_val) < bound(i) + bound(j)) {
        m_targets[std::make_pair(i, j)] = tar_val;
 
        int rel_val;
        if (CGAL::abs(to_lower) < CGAL::abs(to_upper))
          rel_val = ((90 * static_cast<int>(mes_90)) % 180 == 0 ? 0 : 1);
        else
          rel_val = ((90 * static_cast<int>(mes_90 + 1.0)) % 180 == 0 ? 0 : 1);
        
        m_relations[std::make_pair(i, j)] = rel_val;
      } 
      return tar_val;
    }

    /*!
      \brief returns `max_angle`.
    */
    FT bound(const std::size_t) const {
      return m_theta_max;
    }

    /*!
      \brief applies new orientations computed by the QP solver 
      to the initial segments.

      \param solution
      a vector with updated segment orientations.

      \pre `solution.size() > 0`
    */
    void update(const std::vector<FT>& solution) {
      CGAL_precondition(solution.size() > 0);

      Targets_map targets;
      Relations_map relations;
      std::map<std::size_t, Segment_data> segments;
      std::map<FT, Indices> parallel_groups;

      CGAL_precondition(m_targets.size() > 0);
      CGAL_precondition(m_targets.size() == m_relations.size());

      for (const auto& group : m_groups) {
        if (group.size() < 2) continue;

        segments.clear(); targets.clear(); relations.clear();
        build_grouping_data(group, segments, targets, relations);

        parallel_groups.clear();
        if (segments.size() > 0) {
          const std::size_t n = m_input_range.size();

          m_grouping.make_groups(
            m_theta_max, n, segments, solution, 
            parallel_groups, targets, relations);
          rotate_parallel_segments(parallel_groups);
        }
      }
    }

    /// @}

    /// \name Miscellaneous
    /// @{ 

    /*!
      \brief returns indices of parallel segments organized into groups.

      \param groups
      an instance of OutputIterator
    */
    template<typename OutputIterator>
    OutputIterator parallel_groups(OutputIterator groups) {
      for(const auto& parallel_group : m_parallel_groups) {
        const auto& group = parallel_group.second;
        *(groups++) = group;
      }
      return groups;
    }

    /*!
      \brief inserts a group of segments from `input_range`.

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
      
      CGAL_precondition(index_range.size() > 1);
      if (index_range.size() < 2) return;
      
      Indices group;
      group.reserve(index_range.size());
      for (const auto seg_index : index_range)
        group.push_back(seg_index);
      
      m_groups.push_back(group);
      update_segment_data(group);
    }

    /*!
      \brief inserts all input segments from `input_range` as one unique group.

      For more details, 
      see `CGAL::Shape_regularization::Angle_regularization_2::add_group()`.
    */
    void create_unique_group() {
      
      Indices group(m_input_range.size());
      std::iota(group.begin(), group.end(), 0);
      add_group(group);
    }

    /*!
      \brief returns the number of modifed segments.
    */
    std::size_t number_of_modified_segments() const {
      return m_num_modified_segments;
    }

    /// @}

  private:
    Input_range& m_input_range;
    FT m_theta_max;
    const Segment_map m_segment_map;
    std::map<std::size_t, Segment_data> m_segments;
    std::map<Size_pair, FT> m_targets;
    std::map<Size_pair, int> m_relations;
    Grouping m_grouping;
    std::map<FT, Indices> m_parallel_groups;
    std::vector<Indices> m_groups;
    std::size_t m_num_modified_segments;

    void update_segment_data(
      const Indices& group) {
      if (group.size() < 2) return;

      for(const std::size_t seg_index : group) {
        if(m_segments.find(seg_index) != m_segments.end())
          continue;

        const auto& segment = get(m_segment_map, 
          *(m_input_range.begin() + seg_index));
        const Segment_data seg_data(segment, seg_index);
        m_segments.emplace(seg_index, seg_data);
      }
    }

    void build_grouping_data(
      const Indices& group,
      std::map<std::size_t, Segment_data>& segments,
      Targets_map& targets,
      Relations_map& relations) {
      
      for (const std::size_t seg_index : group) {
        CGAL_precondition(m_segments.find(seg_index) != m_segments.end());
        const auto& seg_data = m_segments.at(seg_index);

        segments.emplace(seg_index, seg_data);
        std::size_t tar_index = 0;

        auto rit = m_relations.begin();
        for (const auto& target : m_targets) {
          const std::size_t seg_index_tar_i = target.first.first;
          const std::size_t seg_index_tar_j = target.first.second;
          const FT tar_val = target.second;

          auto& relation = *rit;
          const std::size_t seg_index_rel_i = relation.first.first;
          const std::size_t seg_index_rel_j = relation.first.second;
          const int rel_val = relation.second;

          CGAL_precondition(seg_index_tar_i == seg_index_rel_i);
          CGAL_precondition(seg_index_tar_j == seg_index_rel_j);

          if (seg_index_tar_i == seg_index && seg_index_rel_i == seg_index) {
            targets[std::make_pair(seg_index_tar_i, seg_index_tar_j)] = 
              std::make_pair(tar_val, tar_index);
            relations[std::make_pair(seg_index_rel_i, seg_index_rel_j)] = 
              std::make_pair(rel_val, tar_index);
          }
          ++rit; ++tar_index;
        }
      }
      CGAL_postcondition(targets.size() == relations.size());
    }

    void rotate_parallel_segments(
      const std::map<FT, Indices>& parallel_groups) {
      
      for (const auto& parallel_group : parallel_groups) {
        const FT angle = parallel_group.first;
        const auto& group = parallel_group.second;

        if (m_parallel_groups.find(angle) == m_parallel_groups.end())
          m_parallel_groups[angle] = group;

        // Each group of parallel segments has a normal vector 
        // that we compute with alpha.
        const double angle_rad = 
          CGAL::to_double(angle * static_cast<FT>(CGAL_PI) / FT(180));
        const FT x = static_cast<FT>(std::cos(angle_rad));
        const FT y = static_cast<FT>(std::sin(angle_rad));

        Vector_2 direction = Vector_2(x, y);
        const Vector_2 orth = Vector_2(-direction.y(), direction.x());
        const FT a = orth.x();
        const FT b = orth.y();

        // Rotate segments with precision.
        // Compute equation of the supporting line of the rotated segment.
        for (const std::size_t seg_index : group) {
          CGAL_precondition(m_segments.find(seg_index) != m_segments.end());

          const auto& seg_data = m_segments.at(seg_index);
          const auto& barycenter = seg_data.barycenter;
          const FT c = -a * barycenter.x() - b * barycenter.y();
          set_orientation(seg_index, a, b, c, direction);
        }
      }
    }

    void set_orientation(
      const std::size_t seg_index, 
      const FT a, const FT b, const FT c, 
      Vector_2 direction) {
      
      if (
        direction.y() < FT(0) || (
        direction.y() == FT(0) && direction.x() < FT(0))) 
        direction = -direction;

      FT x1, y1, x2, y2;
      const auto& seg_data = m_segments.at(seg_index);
      const auto& barycenter = seg_data.barycenter;
      const FT length = seg_data.length;

      if (CGAL::abs(direction.x()) > CGAL::abs(direction.y())) { 
        x1 = barycenter.x() - length * direction.x() / FT(2);
        x2 = barycenter.x() + length * direction.x() / FT(2);
        y1 = (-c - a * x1) / b;
        y2 = (-c - a * x2) / b;
      }  else {
        y1 = barycenter.y() - length * direction.y() / FT(2);
        y2 = barycenter.y() + length * direction.y() / FT(2);
        x1 = (-c - b * y1) / a;
        x2 = (-c - b * y2) / a;
      }
      const Point_2 source = Point_2(x1, y1);
      const Point_2 target = Point_2(x2, y2);

      m_input_range[seg_index] = Segment_2(source, target);
      ++m_num_modified_segments;
    } 
  };

} // namespace Segments
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_ANGLE_REGULARIZATION_2_H
