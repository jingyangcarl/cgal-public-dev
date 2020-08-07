#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights/Discrete_harmonic_weight.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using DH = CGAL::Generalized_weights::Discrete_harmonic_weight<Kernel>;

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
  DH dh;
  std::cout << "2D discrete harmonic: " << dh(q2, t2, r2, p2) << std::endl;
  std::cout << "3D discrete harmonic: " << dh(q3, t3, r3, p3) << std::endl;
  std::cout << "-----------------------" << std::endl;

  // Using free functions.
  std::cout << "2D discrete harmonic: " <<
    CGAL::Generalized_weights::discrete_harmonic_weight_2(q2, t2, r2, p2) << std::endl;
  std::cout << "3D discrete harmonic: " <<
    CGAL::Generalized_weights::discrete_harmonic_weight_3(q3, t3, r3, p3) << std::endl;

  return EXIT_SUCCESS;
}
