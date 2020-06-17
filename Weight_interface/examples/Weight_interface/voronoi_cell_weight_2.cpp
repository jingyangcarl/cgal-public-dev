#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Local_averaging_regions_2/Voronoi_cell_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using VC2 = CGAL::Generalized_weights::Voronoi_cell_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2( 0, 0);
  const Point_2 vj2    = Point_2( 0, 1);
  const Point_2 vp2    = Point_2(-1, 0);

  // 3D configuration.
  const Point_3 query3 = Point_3( 0, 0, 1);
  const Point_3 vj3    = Point_3( 0, 1, 1);
  const Point_3 vp3    = Point_3(-1, 0, 1);

  // Compute weights.
  VC2 voronoi;
  std::cout << "2D voronoi cell: " << voronoi(query2, vj2, vp2) << std::endl;
  std::cout << "3D voronoi cell: " << voronoi(query3, vj3, vp3) << std::endl;

  return EXIT_SUCCESS;
}
