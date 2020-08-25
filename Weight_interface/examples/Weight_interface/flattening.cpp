#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/authalic_weights.h>
#include <CGAL/Weight_interface/three_point_family_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

int main() {
  std::cout << std::fixed;

  // 3D configuration.
  const Point_3 q3 = Point_3(3, 1, 1);
  const Point_3 t3 = Point_3(0, 1, 1);
  const Point_3 r3 = Point_3(2, 0, 1);
  const Point_3 p3 = Point_3(7, 1, 1);

  // Choose a type of the weight:
  // e.g. 0 - Wachspress (WP) weight; 1 - mean value (MV);
  const FT wp = FT(0);
  const FT mv = FT(1);

  // Compute WP and MV weights.
  std::cout << "3D wachspress (WP, q3): ";
  std::cout << CGAL::Weights::
    three_point_family_weight(q3, t3, r3, p3, wp) << std::endl;
  std::cout << "3D mean value (MV, q3): ";
  std::cout << CGAL::Weights::
    three_point_family_weight(q3, t3, r3, p3, mv) << std::endl;

  // Converge WP towards MV.
  std::cout << "Converge WP to MV on q3: " << std::endl;
  const FT step = FT(1) / FT(10);
  for (FT x = FT(0); x <= FT(1); x += step) {
    std::cout << "3D x: ";
    std::cout << CGAL::Weights::
      three_point_family_weight(q3, t3, r3, p3, x) << std::endl;
  }

  // Compute WP weights for query3, which is not on the plane [t3, r3, p3].
  Point_3 query3 = Point_3(3, 1, 2);
    std::cout << "3D wachspress (WP, query3): ";
  std::cout << CGAL::Weights::
    three_point_family_weight(query3, t3, r3, p3, wp) << std::endl;

  // Converge query3 towards q3 that is we flatten the configuration.
  // We also compare the result with the authalic weight.
  std::cout << "Converge query3 to q3: " << std::endl;
  for (FT x = FT(0); x <= FT(1); x += step) {
    std::cout << "3D wachspress/authalic: ";
    query3 = Point_3(3, 1, FT(2) - x);
    std::cout << CGAL::Weights::
      three_point_family_weight(query3, t3, r3, p3, wp) << "/";
    std::cout << CGAL::Weights::
      authalic_weight(query3, t3, r3, p3) << std::endl;
  }

  return EXIT_SUCCESS;
}
