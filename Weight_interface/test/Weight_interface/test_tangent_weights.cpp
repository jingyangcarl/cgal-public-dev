#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights/utils.h>
#include <CGAL/Weight_interface/Generalized_weights/tangent_weights.h>

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
  std::cout << "2D tangent: " <<
    CGAL::Generalized_weights::tangent_weight_2(q2, t2, r2, p2) << std::endl;
  std::cout << "3D tangent: " <<
    CGAL::Generalized_weights::tangent_weight_3(q3, t3, r3, p3) << std::endl;
  std::cout << "-------------" << std::endl;

  // Construct a 2D weight.
  const auto w2 =
    CGAL::Generalized_weights::half_tangent_weight(
      CGAL::Generalized_weights::utils::distance_2(r2, q2),
      CGAL::Generalized_weights::utils::distance_2(t2, q2),
      CGAL::Generalized_weights::utils::area_2(r2, q2, t2),
      CGAL::Generalized_weights::utils::scalar_product_2(r2, q2, t2)) +
    CGAL::Generalized_weights::half_tangent_weight(
      CGAL::Generalized_weights::utils::distance_2(r2, q2),
      CGAL::Generalized_weights::utils::distance_2(p2, q2),
      CGAL::Generalized_weights::utils::area_2(p2, q2, r2),
      CGAL::Generalized_weights::utils::scalar_product_2(p2, q2, r2));
  std::cout << "2D tangent: " << w2 << std::endl;

  // Construct a 3D weight.
  const auto w3 =
    CGAL::Generalized_weights::half_tangent_weight(
      CGAL::Generalized_weights::tangent_half_angle(
        CGAL::Generalized_weights::utils::distance_3(r3, q3),
        CGAL::Generalized_weights::utils::distance_3(t3, q3),
        CGAL::Generalized_weights::utils::area_3(r3, q3, t3),
        CGAL::Generalized_weights::utils::scalar_product_3(r3, q3, t3)),
      CGAL::Generalized_weights::utils::distance_3(r3, q3)) +
    CGAL::Generalized_weights::half_tangent_weight(
      CGAL::Generalized_weights::tangent_half_angle(
        CGAL::Generalized_weights::utils::distance_3(r3, q3),
        CGAL::Generalized_weights::utils::distance_3(p3, q3),
        CGAL::Generalized_weights::utils::area_3(p3, q3, r3),
        CGAL::Generalized_weights::utils::scalar_product_3(p3, q3, r3)),
      CGAL::Generalized_weights::utils::distance_3(r3, q3));
  std::cout << "3D tangent: " << w3 << std::endl;

  return EXIT_SUCCESS;
}
