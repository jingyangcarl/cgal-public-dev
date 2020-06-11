#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Uniform_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using UN2 = CGAL::Generalized_weights::Uniform_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2(+FT(0), FT(0));
  const Point_2 vm2    = Point_2(+FT(1), FT(0));
  const Point_2 vi2    = Point_2(+FT(0), FT(1));
  const Point_2 vp2    = Point_2(-FT(1), FT(0));

  // 3D configuration.
  const Point_3 query3 = Point_3(+FT(0), FT(0), FT(1));
  const Point_3 vm3    = Point_3(+FT(1), FT(0), FT(1));
  const Point_3 vi3    = Point_3(+FT(0), FT(1), FT(1));
  const Point_3 vp3    = Point_3(-FT(1), FT(0), FT(1));

  // Polygon mesh configuration.
  // add polygon_mesh, vdi, and vcj

  // Compute weights.
  UN2 uniform;
  std::cout << "2D: " << uniform(query2, vm2, vi2, vp2) << std::endl;
  std::cout << "3D: " << uniform(query3, vm3, vi3, vp3) << std::endl;
  // std::cout << "PM: " << uniform(polygon_mesh, vdi, vcj) << std::endl;

  return EXIT_SUCCESS;
}
