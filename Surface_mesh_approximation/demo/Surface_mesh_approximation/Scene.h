#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/internal/Surface_mesh_approximation/VSA_segmentation.h>
#include "types.h"

class Scene
{
public:
  // types
  typedef CGAL::Bbox_3 Bbox;
  typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type PointPropertyMap;
  typedef CGAL::internal::VSA_segmentation<Polyhedron, Kernel, PointPropertyMap> VSA;
  typedef VSA::Anchor Anchor;

public:
  Scene();
  ~Scene();

  void update_bbox();
  Bbox bbox() { return m_bbox; }

  // file menu
  int open(QString filename);

  // algorithms
  void VSA_segmentation(const std::size_t num_proxies, const std::size_t num_iterations);
  void VSA_incremental(const std::size_t num_proxies, const std::size_t num_iterations);

  // toggle view options
  void toggle_view_wireframe() {
    m_view_wireframe = !m_view_wireframe;
  }

  void toggle_view_seg_boundary() {
    m_view_seg_boundary = !m_view_seg_boundary;
  }

  void toggle_view_anchors() {
    m_view_anchors = !m_view_anchors;
  }

  void toggle_view_approximation() {
    m_view_approximation = !m_view_approximation;
  }

  void draw();

private:
  Vector normalize(const Vector& v) {
    return v / std::sqrt(v * v);
  }

  // rendering
  void render_polyhedron();
  void render_wireframe();
  void render_segment_boundary();
  void render_anchors(const double offset = 0);
  void render_borders(const double offset = 0);
  void render_approximation(const double);

private:
  // member data
  Bbox m_bbox;
  Polyhedron *m_pPolyhedron;

  typedef std::map<Polyhedron::Facet_const_handle, std::size_t> FacetIdMap;
  typedef boost::associative_property_map<FacetIdMap> FacetIdPmap;
  FacetIdMap m_fidx_map;
  FacetIdPmap m_fidx_pmap; // property-map for segment-idx

  std::vector<Anchor> m_anchors; // the anchors
  std::vector<std::vector<std::size_t> > m_bdrs; // anchor borders
  std::vector<int> m_tris;

  std::size_t m_px_num;

  // view options
  bool m_view_wireframe;
  bool m_view_seg_boundary;
  bool m_view_anchors;
  bool m_view_approximation;
}; // end class Scene

#endif // SCENE_H