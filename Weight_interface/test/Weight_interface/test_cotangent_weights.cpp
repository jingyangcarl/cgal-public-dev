#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights/utils.h>
#include <CGAL/Weight_interface/Generalized_weights/cotangent_weights.h>

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
  std::cout << "2D cotangent: " <<
    CGAL::Generalized_weights::cotangent_weight(q2, t2, r2, p2) << std::endl;
  std::cout << "3D cotangent: " <<
    CGAL::Generalized_weights::cotangent_weight(q3, t3, r3, p3) << std::endl;
  std::cout << "---------------" << std::endl;

  // Construct a 2D weight.
  const auto w2 =
    CGAL::Generalized_weights::half_cotangent_weight(
      CGAL::Generalized_weights::utils::cotangent(q2, t2, r2)) +
    CGAL::Generalized_weights::half_cotangent_weight(
      CGAL::Generalized_weights::utils::cotangent(r2, p2, q2));
  std::cout << "2D cotangent: " << w2 << std::endl;

  // Construct a 3D weight.
  const auto w3 =
    CGAL::Generalized_weights::half_cotangent_weight(
      CGAL::Generalized_weights::utils::cotangent(q3, t3, r3)) +
    CGAL::Generalized_weights::half_cotangent_weight(
      CGAL::Generalized_weights::utils::cotangent(r3, p3, q3));
  std::cout << "3D cotangent: " << w3 << std::endl;

  return EXIT_SUCCESS;
}
