#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights/utils.h>
#include <CGAL/Weight_interface/Generalized_weights/cotangent_weights.h>
#include <CGAL/Weight_interface/Weighting_regions/mixed_voronoi_region_weights.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Surface_mesh.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Solver_traits = CGAL::Eigen_solver_traits<
  Eigen::SparseLU<CGAL::Eigen_sparse_matrix<FT>::EigenType> >;
using Matrix = Solver_traits::Matrix;

using Mesh = CGAL::Surface_mesh<Point_3>;
using VD   = boost::graph_traits<Mesh>::vertex_descriptor;
using HD   = boost::graph_traits<Mesh>::halfedge_descriptor;

using Indices = std::map<VD, std::size_t>;
using VIMap   = boost::associative_property_map<Indices>;

struct Weight_wrapper {

  template<typename PointMap>
  FT w_ij(
    const Mesh& mesh, const HD he, const PointMap pmap) const {

    const auto v0 = target(he, mesh);
    const auto v1 = source(he, mesh);

    const auto& p0 = get(pmap, v0);
    const auto& p1 = get(pmap, v1);

    if (is_border_edge(he, mesh)) {
      const auto he_cw = opposite(next(he, mesh), mesh);
      auto v2 = source(he_cw, mesh);

      if (is_border_edge(he_cw, mesh)) {
        const auto he_ccw = prev(opposite(he, mesh), mesh);
        v2 = source(he_ccw, mesh);

        const auto& p2 = get(pmap, v2);
        return CGAL::Generalized_weights::utils::cotangent(p1, p2, p0);
      } else {
        const auto& p2 = get(pmap, v2);
        return CGAL::Generalized_weights::utils::cotangent(p0, p2, p1);
      }
    }

    const auto he_cw = opposite(next(he, mesh), mesh);
    const auto v2 = source(he_cw, mesh);
    const auto he_ccw = prev(opposite(he, mesh), mesh);
    const auto v3 = source(he_ccw, mesh);

    const auto& p2 = get(pmap, v2);
    const auto& p3 = get(pmap, v3);
    return CGAL::Generalized_weights::
      cotangent_weight(p0, p2, p1, p3) / FT(2);
  }

  template<typename PointMap>
  FT w_i(
    const Mesh& mesh, const VD v_i, const PointMap pmap) const {

    FT A_i = FT(0);
    const auto v0 = v_i;
    const auto init = halfedge(v_i, mesh);
    for (const auto he : halfedges_around_target(init, mesh)) {
      assert(v0 == target(he, mesh));
      if (is_border(he, mesh)) continue;

      const auto v1 = source(he, mesh);
      const auto v2 = target(next(he, mesh), mesh);

      const auto& p0 = get(pmap, v0);
      const auto& p1 = get(pmap, v1);
      const auto& p2 = get(pmap, v2);

      A_i += CGAL::Generalized_weights::
        mixed_voronoi_area(p1, p0, p2);
    }
    assert(A_i != FT(0));
    return FT(1) / (FT(2) * A_i);
  }
};

void set_laplacian_matrix(
  const Mesh& mesh, const VIMap& vimap, Matrix& L) {

  const Weight_wrapper wrapper;
  const auto pmap = get(CGAL::vertex_point, mesh);

  // Precompute Voronoi areas.
  std::map<std::size_t, FT> w_i;
  for (const auto v_i : vertices(mesh))
    w_i[get(vimap, v_i)] = wrapper.w_i(mesh, v_i, pmap);

  // Fill the matrix.
  for (const auto he : halfedges(mesh)) {
    const auto vi = source(he, mesh);
    const auto vj = target(he, mesh);

    const std::size_t i = get(vimap, vi);
    const std::size_t j = get(vimap, vj);
    if (i > j) continue;

    const FT w_ij = w_i[j] * wrapper.w_ij(mesh, he, pmap);
    L.set_coef(i, j,  w_ij, true);
    L.set_coef(j, i,  w_ij, true);
    L.add_coef(i, i, -w_ij);
    L.add_coef(j, j, -w_ij);
  }
}

int main() {

  // Create mesh.
  Mesh mesh;
  const auto v0 = mesh.add_vertex(Point_3(0, 2, 0));
  const auto v1 = mesh.add_vertex(Point_3(2, 2, 0));
  const auto v2 = mesh.add_vertex(Point_3(0, 0, 0));
  const auto v3 = mesh.add_vertex(Point_3(2, 0, 0));
  const auto v4 = mesh.add_vertex(Point_3(1, 1, 1));
  mesh.add_face(v0, v2, v4);
  mesh.add_face(v2, v3, v4);
  mesh.add_face(v3, v1, v4);
  mesh.add_face(v1, v0, v4);
  mesh.add_face(v2, v3, v1);
  mesh.add_face(v1, v0, v2);
  assert(CGAL::is_triangle_mesh(mesh));

  // Create vertex to index map.
  Indices indices;
  VIMap vimap(indices);
  std::size_t n = 0;
  for (const auto vd : vertices(mesh))
    put(vimap, vd, n++);
  assert(n == 5);

  // Set discretized Laplacian.
  Matrix L(n, n);
  set_laplacian_matrix(mesh, vimap, L);
  std::cout << std::fixed; std::cout << std::showpos;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j)
      std::cout << L.get_coef(i, j) << " ";
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
