#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

// #define UF_WEIGHT
// #define ID_WEIGHT
// #define SP_WEIGHT
// #define TP_WEIGHT
// #define WP_WEIGHT
#define MV_WEIGHT
// #define DH_WEIGHT

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
  std::cout << "2D/3D weight: ";
  #if defined(UF_WEIGHT)
    std::cout << CGAL::Generalized_weights::uniform_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::uniform_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(ID_WEIGHT)
    std::cout << CGAL::Generalized_weights::inverse_distance_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::inverse_distance_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(SP_WEIGHT)
    std::cout << CGAL::Generalized_weights::shepard_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::shepard_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(TP_WEIGHT)
    std::cout << CGAL::Generalized_weights::three_point_family_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::three_point_family_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(WP_WEIGHT)
    std::cout << CGAL::Generalized_weights::wachspress_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::wachspress_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(MV_WEIGHT)
    std::cout << CGAL::Generalized_weights::mean_value_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::mean_value_weight_3(q3, t3, r3, p3) << std::endl;
  #elif defined(DH_WEIGHT)
    std::cout << CGAL::Generalized_weights::discrete_harmonic_weight_2(q2, t2, r2, p2) << "/";
    std::cout << CGAL::Generalized_weights::discrete_harmonic_weight_3(q3, t3, r3, p3) << std::endl;
  #endif

  return EXIT_SUCCESS;
}
