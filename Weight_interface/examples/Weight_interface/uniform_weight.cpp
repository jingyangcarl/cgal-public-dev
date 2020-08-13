#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights/uniform_weight.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

int main() {

  // 2D configuration.
  const Point_2 q2 = Point_2( 0,  0);
  const Point_2 t2 = Point_2(-1,  0);
  const Point_2 r2 = Point_2( 0, -1);
  const Point_2 p2 = Point_2( 1,  0);

  // 3D configuration.
  const Point_3 q3 = Point_3( 0,  0, 1);
  const Point_3 t3 = Point_3(-1,  0, 1);
  const Point_3 r3 = Point_3( 0, -1, 1);
  const Point_3 p3 = Point_3( 1,  0, 1);

  // Compute weights.
  std::cout << "2D uniform: " <<
    CGAL::Generalized_weights::uniform_weight_2(q2, t2, r2, p2) << std::endl;
  std::cout << "3D uniform: " <<
    CGAL::Generalized_weights::uniform_weight_3(q3, t3, r3, p3) << std::endl;
  std::cout << "-------------" << std::endl;

  // Overloads.
  std::cout << "2D uniform: " <<
    CGAL::Generalized_weights::uniform_weight_2() << std::endl;
  std::cout << "3D uniform: " <<
    CGAL::Generalized_weights::uniform_weight_3() << std::endl;

  return EXIT_SUCCESS;
}
