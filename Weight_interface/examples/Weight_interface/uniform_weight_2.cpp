#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Uniform_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using UW2 = CGAL::Generalized_weights::Uniform_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2( 0, 0);
  const Point_2 vm2    = Point_2( 1, 0);
  const Point_2 vj2    = Point_2( 0, 1);
  const Point_2 vp2    = Point_2(-1, 0);

  // 3D configuration.
  const Point_3 query3 = Point_3( 0, 0, 1);
  const Point_3 vm3    = Point_3( 1, 0, 1);
  const Point_3 vj3    = Point_3( 0, 1, 1);
  const Point_3 vp3    = Point_3(-1, 0, 1);

  // Compute weights.
  UW2 uniform;
  std::cout << "2D uniform: " << uniform(query2, vm2, vj2, vp2) << std::endl;
  std::cout << "3D uniform: " << uniform(query3, vm3, vj3, vp3) << std::endl;

  // 2D half configuration.
  const Point_2 queryh2 = Point_2( 0, 0);
  const Point_2 vjh2    = Point_2( 0, 1);
  const Point_2 vph2    = Point_2(-1, 0);

  // 3D half configuration.
  const Point_3 queryh3 = Point_3( 0, 0, 1);
  const Point_3 vjh3    = Point_3( 0, 1, 1);
  const Point_3 vph3    = Point_3(-1, 0, 1);

  // Compute half weights.
  std::cout << "2D uniform half: " << uniform(queryh2, vjh2, vph2) << std::endl;
  std::cout << "3D uniform half: " << uniform(queryh3, vjh3, vph3) << std::endl;

  return EXIT_SUCCESS;
}
