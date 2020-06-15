#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Authalic_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using AW2 = CGAL::Generalized_weights::Authalic_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2(+FT(0), FT(0));
  const Point_2 vm2    = Point_2(+FT(1), FT(0));
  const Point_2 vj2    = Point_2(+FT(0), FT(1));
  const Point_2 vp2    = Point_2(-FT(1), FT(0));

  // 3D configuration.
  const Point_3 query3 = Point_3(+FT(0), FT(0), FT(1));
  const Point_3 vm3    = Point_3(+FT(1), FT(0), FT(1));
  const Point_3 vj3    = Point_3(+FT(0), FT(1), FT(1));
  const Point_3 vp3    = Point_3(-FT(1), FT(0), FT(1));

  // Compute weights.
  AW2 authalic;
  std::cout << "2D authalic: " << authalic(query2, vm2, vj2, vp2) << std::endl;
  std::cout << "3D authalic: " << authalic(query3, vm3, vj3, vp3) << std::endl;

  return EXIT_SUCCESS;
}
