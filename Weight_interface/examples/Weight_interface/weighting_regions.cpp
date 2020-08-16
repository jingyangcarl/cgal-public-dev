#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Weighting_regions.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

// #define UF_WEIGHT
// #define TR_WEIGHT
// #define BC_WEIGHT
#define VN_WEIGHT
// #define MV_WEIGHT

int main() {

  // 2D configuration.
  const Point_2 p2 = Point_2(-1, 0);
  const Point_2 q2 = Point_2( 0, 0);
  const Point_2 r2 = Point_2( 0, 1);

  // 3D configuration.
  const Point_3 p3 = Point_3(-1, 0, 1);
  const Point_3 q3 = Point_3( 0, 0, 1);
  const Point_3 r3 = Point_3( 0, 1, 1);

  // Compute weights.
  std::cout << "2D/3D weight: ";
  #if defined(UF_WEIGHT)
    std::cout << CGAL::Generalized_weights::uniform_area_2(p2, q2, r2) << "/";
    std::cout << CGAL::Generalized_weights::uniform_area_3(p3, q3, r3) << std::endl;
  #elif defined(TR_WEIGHT)
    std::cout << CGAL::Generalized_weights::triangle_area_2(p2, q2, r2) << "/";
    std::cout << CGAL::Generalized_weights::triangle_area_3(p3, q3, r3) << std::endl;
  #elif defined(BC_WEIGHT)
    std::cout << CGAL::Generalized_weights::barycentric_area_2(p2, q2, r2) << "/";
    std::cout << CGAL::Generalized_weights::barycentric_area_3(p3, q3, r3) << std::endl;
  #elif defined(VN_WEIGHT)
    std::cout << CGAL::Generalized_weights::voronoi_area_2(p2, q2, r2) << "/";
    std::cout << CGAL::Generalized_weights::voronoi_area_3(p3, q3, r3) << std::endl;
  #elif defined(MV_WEIGHT)
    std::cout << CGAL::Generalized_weights::mixed_voronoi_area_2(p2, q2, r2) << "/";
    std::cout << CGAL::Generalized_weights::mixed_voronoi_area_3(p3, q3, r3) << std::endl;
  #endif

  return EXIT_SUCCESS;
}
