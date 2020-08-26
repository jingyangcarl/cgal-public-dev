#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

int main() {

  // 2D configuration.
  const Point_2 t2 = Point_2(-1,  0); // p0
  const Point_2 r2 = Point_2( 0, -1); // p1
  const Point_2 p2 = Point_2( 1,  0); // p2
  const Point_2 q2 = Point_2( 0,  0); // query

  // 3D configuration.
  const Point_3 t3 = Point_3(-1,  0, 1); // p0
  const Point_3 r3 = Point_3( 0, -1, 1); // p1
  const Point_3 p3 = Point_3( 1,  0, 1); // p2
  const Point_3 q3 = Point_3( 0,  0, 1); // query

  std::cout << "2D/3D tangent weight: ";
  std::cout << CGAL::Weights::tangent_weight(t2, r2, p2, q2) << "/";
  std::cout << CGAL::Weights::tangent_weight(t3, r3, p3, q3) << std::endl;

  std::cout << "2D/3D shepard weight: ";
  std::cout << CGAL::Weights::shepard_weight(t2, r2, p2, q2, 2.0) << "/";
  std::cout << CGAL::Weights::shepard_weight(t3, r3, p3, q3, 2.0) << std::endl;

  std::cout << "2D/3D barycentric area: ";
  std::cout << CGAL::Weights::barycentric_area(p2, q2, r2) << "/";
  std::cout << CGAL::Weights::barycentric_area(p3, q3, r3) << std::endl;

  return EXIT_SUCCESS;
}
