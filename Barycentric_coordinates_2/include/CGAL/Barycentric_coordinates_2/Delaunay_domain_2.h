// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Dmitry Anisimov, David Bommes, Kai Hormann, Pierre Alliez
//

#ifndef CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
#define CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// CGAL includes.
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>
#include <CGAL/Barycentric_coordinates_2/barycentric_enum_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefHarmonic

    \brief 2D Delaunay domain restricted to a simple polygon.

    This class implements a discretized domain restricted to a simple polygon.
    The interior part of the input polygon is triangulated and refined with respect
    to the user-defined shape size parameter. The final triangulation is a constrained
    Delaunay triangulation, where the constraints are the polygon edges.

    Internally, the package \ref PkgMesh2 is used. See it for more details.

    \tparam Polygon
    is a model of `ConstRange`.

    \tparam GeomTraits
    is a model of `BarycentricTraits_2`.

    \tparam VertexMap
    is a `ReadablePropertyMap` whose key type is `Polygon::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \cgalModels `DiscretizedDomain_2`
  */
  template<
  typename Polygon,
  typename GeomTraits,
  typename VertexMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Delaunay_domain_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Polygon_ = Polygon;
    using GT = GeomTraits;
    using Vertex_map = VertexMap;

    struct VI {
      bool is_on_boundary = false;
      std::size_t index = std::size_t(-1);
      std::vector<std::size_t> neighbors;
    };

    using FB  = CGAL::Delaunay_mesh_face_base_2<GeomTraits>;
    using VB  = CGAL::Triangulation_vertex_base_with_info_2<VI, GeomTraits>;
    using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<GeomTraits, TDS, TAG>;

    using Vertex_handle = typename CDT::Vertex_handle;

    using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
    using Mesher   = CGAL::Delaunay_mesher_2<CDT, Criteria>;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param polygon
      An instance of `Polygon` with the vertices of a simple polygon.

      \param traits
      An instance of `GeomTraits`. The default initialization is provided.

      \param vertex_map
      An instance of `VertexMap` that maps a vertex from `polygon`
      to `Point_2`. The default is the identity property map.

      \pre `polygon.size() >= 3`
      \pre `polygon is simple`
    */
    Delaunay_domain_2(
      const Polygon& polygon,
      const GeomTraits traits = GeomTraits(),
      const VertexMap vertex_map = VertexMap()) :
    m_traits(traits) {

      m_polygon.clear();
      m_polygon.reserve(polygon.size());
      for (const auto& item : polygon)
        m_polygon.push_back(get(vertex_map, item));
      CGAL_precondition(m_polygon.size() >= 3);

      CGAL_precondition(
        CGAL::is_simple_2(m_polygon.begin(), m_polygon.end(), m_traits));
    }

    /*!
      \brief creates a refined Delaunay triangulation restricted to the input polygon.

      After the construction is completed, the first n vertices are the polygon vertices,
      the next m vertices are the vertices generated along the polygon boundary,
      the last k vertices are the vertices generated in the interior part of the polygon.

      \param shape_size
      A shape size bound. See `Delaunay_mesh_size_criteria_2` for more details.

      \param list_of_seeds
      Contains seed points indicating, which parts of the input polygon
      should be partitioned and subdivided.
    */
    void create(
      const FT shape_size,
      const std::list<Point_2>& list_of_seeds) {

      create_triangulation();
      refine_triangulation(
        shape_size, list_of_seeds);
      check_boundaries();
      create_neighbors();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes barycenters of all generated triangles.

      \param barycenters
      An `std::vector` that stores the computed barycenters.
    */
    void barycenters(
      std::vector<Point_2>& barycenters) const {

      const std::size_t num_faces = get_number_of_faces();
      if (num_faces == 0) return;

      barycenters.clear();
      barycenters.reserve(num_faces);

      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        const Point_2 b = CGAL::barycenter(
        fh->vertex(0)->point(), FT(1),
        fh->vertex(1)->point(), FT(1),
        fh->vertex(2)->point(), FT(1));
        barycenters.push_back(b);
      }
      CGAL_assertion(barycenters.size() == num_faces);
    }

    /*!
      \brief returns the number of triangulation vertices.

      This function implements `DiscretizedDomain_2::number_of_vertices()`.
    */
    const std::size_t number_of_vertices() const {

      CGAL_assertion(
        m_vhs.size() == m_cdt.number_of_vertices());
      return m_vhs.size();
    }

    /*!
      \brief returns a const reference to the triangulation vertex with
      the index `query_index`.

      This function implements `DiscretizedDomain_2::vertex()`.

      \param query_index
      An index of the requested vertex.

      \pre `query_index >= 0 && query_index < number_of_vertices()`
    */
    const Point_2& vertex(
      const std::size_t query_index) const {

      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());
      return m_vhs[query_index]->point();
    }

    /*!
      \brief controls if the triangulation vertex with the index `query_index`
      is on the polygon boundary.

      This function implements `DiscretizedDomain_2::is_on_boundary()`.

      \param query_index
      An index of the query vertex.

      \pre `query_index >= 0 && query_index < number_of_vertices()`
    */
    const bool is_on_boundary(
      const std::size_t query_index) const {

      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());
      return m_vhs[query_index]->info().is_on_boundary;
    }

    /*!
      \brief returns the one-ring neighborhood of the triangulation vertex
      with the index `query_index`.

      This function implements `DiscretizedDomain_2::operator()()`.

      \param query_index
      An index of the query vertex.

      \param neighbors
      Stores indices of the vertices, which from the one-ring neighborhood
      of the query vertex.

      \pre `query_index >= 0 && query_index < number_of_vertices()`
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(
        query_index >= 0 && query_index < number_of_vertices());
      const auto vh = m_vhs[query_index];
      neighbors = vh->info().neighbors;
    }

    /*!
      \brief locates a triangle that contains a given query point.

      If `triangle` is empty, the query point does not belong to the domain.

      This function implements `DiscretizedDomain_2::locate()`.

      \param query
      A query point.

      \param triangle
      Stores indices of the vertices, which form a triangle, that contains
      the query point.
    */
    void locate(
      const Point_2& query,
      std::vector<std::size_t>& triangle) const {

      triangle.clear();
      const auto fh = m_cdt.locate(query);
      if (fh->is_in_domain()) {
        for (std::size_t i = 0; i < 3; ++i)
          triangle.push_back(fh->vertex(i)->info().index);
      }
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_vhs.clear();
      m_cdt.clear();
    }

    /*!
      \brief releases all memory that is used internally.
    */
    void release_memory() {
      clear();
      m_vhs.shrink_to_fit();
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void export_points_2(
      const std::vector<Point_2>& points,
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_vertices = points.size();
      add_ply_header_points(num_vertices, out);

      for (const auto& point : points)
        out << point << " 0 0 0 0" << std::endl;
      save(out, file_path + ".ply");
    }

    template<typename Point_3>
    void export_points_3(
      const std::vector<Point_3>& points,
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_vertices = points.size();
      add_ply_header_points(num_vertices, out);

      for (const auto& point : points)
        out << point << " 0 0 0" << std::endl;
      save(out, file_path + ".ply");
    }

    void export_triangulation(
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_faces = get_number_of_faces();
      add_ply_header_mesh(num_faces * 3, num_faces, out);

      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        out << p0.x() << " " << p0.y() << " 0" << std::endl;
        out << p1.x() << " " << p1.y() << " 0" << std::endl;
        out << p2.x() << " " << p2.y() << " 0" << std::endl;
      }

      std::size_t i = 0;
      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        out << 3 << " "
        << i + 0 << " " << i + 1 << " " << i + 2 << " "
        << "0 0 0" << std::endl;
        i += 3;
      }
      save(out, file_path + ".ply");
    }
    /// \endcond

  private:

    // Fields.
    const GeomTraits m_traits;

    CDT m_cdt;
    std::vector<Vertex_handle> m_vhs;

    std::vector<Point_2> m_polygon;

    void create_triangulation() {

      m_cdt.clear(); m_vhs.clear();
      m_vhs.reserve(m_polygon.size());
      for (const auto& vertex : m_polygon)
        m_vhs.push_back(m_cdt.insert(vertex));

      CGAL_assertion(m_vhs.size() == m_polygon.size());
      CGAL_assertion(m_vhs.size() == m_cdt.number_of_vertices());

      for (std::size_t i = 0; i < m_vhs.size(); ++i) {
        const std::size_t ip = (i + 1) % m_vhs.size();
        if (m_vhs[i] != m_vhs[ip])
          m_cdt.insert_constraint(m_vhs[i], m_vhs[ip]);
      }
    }

    void refine_triangulation(
      const FT shape_size,
      const std::list<Point_2>& list_of_seeds) {

      Mesher mesher(m_cdt);
      mesher.set_seeds(list_of_seeds.begin(), list_of_seeds.end(), true);
      mesher.set_criteria(Criteria(shape_size, shape_size));
      mesher.refine_mesh();

      m_vhs.clear();
      m_vhs.reserve(m_cdt.number_of_vertices());

      std::size_t count = 0;
      for (auto vh = m_cdt.finite_vertices_begin();
      vh != m_cdt.finite_vertices_end(); ++vh, ++count) {
        vh->info().index = count;
        m_vhs.push_back(vh);
      }
      CGAL_assertion(m_vhs.size() == m_cdt.number_of_vertices());
    }

    void check_boundaries() {

      for (std::size_t i = 0; i < m_vhs.size(); ++i) {
        const auto vh = m_vhs[i];
        vh->info().is_on_boundary = false;

        auto face = m_cdt.incident_faces(vh);
        CGAL_assertion(!face.is_empty());
        const auto end = face;

        do {
          if (!face->is_in_domain()) {
            vh->info().is_on_boundary = true; break;
          } ++face;
        } while (face != end);
      }
    }

    void create_neighbors() {

      for (std::size_t i = 0; i < m_vhs.size(); ++i) {
        const auto vh = m_vhs[i];

        auto edge = m_cdt.incident_edges(vh);
        CGAL_assertion(!edge.is_empty());
        auto end = edge;

        // Move the pointer to the most left boundary face.
        if (vh->info().is_on_boundary) {
          do {
            const auto fh = edge->first;
            if (!fh->is_in_domain()) {
              end = edge; break;
            } ++edge;
          } while (edge != end);
          do {
            const auto fh = edge->first;
            if (fh->is_in_domain()) {
              end = edge; break;
            } ++edge;
          } while (edge != end);
        }

        // We then traverse edges in the counterclockwise order.
        vh->info().neighbors.clear();
        do {
          const auto fh = edge->first;
          const auto nh = fh->vertex(edge->second);

          if (fh->is_in_domain())
            vh->info().neighbors.push_back(nh->info().index);
          else {
            vh->info().neighbors.push_back(nh->info().index);
            break;
          } ++edge;
        } while (edge != end);
        CGAL_assertion(
          vh->info().neighbors.size() > 0);
      }
    }

    const std::size_t get_number_of_faces() const {

      std::size_t num_faces = 0;
      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;
        ++num_faces;
      }
      return num_faces;
    }

    void add_ply_header_points(
      const std::size_t num_vertices,
      std::stringstream& out) const {

      out <<
			"ply" << std::endl <<
			"format ascii 1.0" << std::endl <<
			"element vertex " << num_vertices << std::endl <<
			"property double x" << std::endl <<
			"property double y" << std::endl <<
			"property double z" << std::endl <<
			"property uchar red"   << std::endl <<
			"property uchar green" << std::endl <<
			"property uchar blue"  << std::endl <<
			"end_header" << std::endl;
    }

    void add_ply_header_mesh(
      const std::size_t num_vertices,
      const std::size_t num_faces,
      std::stringstream& out) const {

      out <<
			"ply" << std::endl <<
			"format ascii 1.0" << std::endl <<
			"element vertex " << num_vertices << std::endl <<
			"property double x" << std::endl <<
			"property double y" << std::endl <<
			"property double z" << std::endl <<
			"element face " << num_faces << std::endl <<
			"property list uchar int vertex_indices" << std::endl <<
			"property uchar red"   << std::endl <<
			"property uchar green" << std::endl <<
			"property uchar blue"  << std::endl <<
			"end_header" << std::endl;
    }

    void save(
      const std::stringstream& out,
      const std::string file_path) const {

      std::ofstream file(
        file_path.c_str(), std::ios_base::out);
      if (!file) {
        std::cerr << std::endl <<
          "ERROR: Error saving file " << file_path
        << "!" << std::endl << std::endl;
      }
      file << out.str();
      file.close();
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
