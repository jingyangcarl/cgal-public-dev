#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Local_averaging_regions_2/Triangle_cell_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using TC2 = CGAL::Generalized_weights::Triangle_cell_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2(+FT(0), FT(0));
  const Point_2 vj2    = Point_2(+FT(0), FT(1));
  const Point_2 vp2    = Point_2(-FT(1), FT(0));

  // 3D configuration.
  const Point_3 query3 = Point_3(+FT(0), FT(0), FT(1));
  const Point_3 vj3    = Point_3(+FT(0), FT(1), FT(1));
  const Point_3 vp3    = Point_3(-FT(1), FT(0), FT(1));

  // Compute weights.
  TC2 triangle;
  std::cout << "2D triangle cell: " << triangle(query2, vj2, vp2) << std::endl;
  std::cout << "3D triangle cell: " << triangle(query3, vj3, vp3) << std::endl;

  return EXIT_SUCCESS;
}
