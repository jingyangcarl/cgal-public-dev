#include <CGAL/Simple_cartesian.h>
#include <CGAL/Weight_interface/Generalized_weights_2/Wachspress_weight_2.h>

// Typedefs.
using Kernel  = CGAL::Simple_cartesian<double>;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Point_3 = typename Kernel::Point_3;

using WP2 = CGAL::Generalized_weights::Wachspress_weight_2<Kernel>;

int main() {

  // 2D configuration.
  const Point_2 query2 = Point_2(+FT(0), FT(0));
  const Point_2 vm2    = Point_2(+FT(1), FT(0));
  const Point_2 vj2    = Point_2(+FT(0), FT(1));
  const Point_2 vp2    = Point_2(-FT(1), FT(0));
  // Am = 0.5, Aj = 0.5, C = 1, B = 0.

  // 3D configuration 1. All points are coplanar on a horizontal plane.
  const Point_3 query3 = Point_3(+FT(0), FT(0), FT(1));
  const Point_3 vm3    = Point_3(+FT(1), FT(0), FT(1));
  const Point_3 vj3    = Point_3(+FT(0), FT(1), FT(1));
  const Point_3 vp3    = Point_3(-FT(1), FT(0), FT(1));
  // Am = 0.5, Aj = 0.5, C = 1, B = 0, 2D : 4.

  // 3D configuration 2. All points are coplanar on an arbitrary plane.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.812058795517);
  // const Point_3 vp3    = Point_3(+4.54574179746300, +7.9148892979530, +0.000000000000);
  // Am = 17.701, Aj = 26.574, C = 13.804, B = 30.471, 3D: 0.0293453.

  // 3D configuration 3. The first three points are coplanar whereas vp3 is slightly offsetted.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.0000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.0000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.8120587955170);
  // const Point_3 vp3    = Point_3(+4.22658807629900, +8.2522664687610, -0.2914612626125);

  // My area Am should coincide with the area below, but my areas Aj, C, and B should be
  // closer to the second configuration than areas below.
  // Am = 17.701, Aj = 26.641, C = 13.896, B = 30.524, 3D: 0.0293359.

  // 3D configuration 3. The first three points are coplanar whereas vp3 is even more slightly offsetted.
  // This configuration should converge to the second one comparing to the previous one.
  // const Point_3 query3 = Point_3(-3.25022842347400, +3.9956322466210, +4.0000000000000);
  // const Point_3 vm3    = Point_3(-4.82908178925400, -0.9535525880045, +0.0000000000000);
  // const Point_3 vj3    = Point_3(-0.09833607144587, +1.9562179601920, -1.8120587955170);
  // const Point_3 vp3    = Point_3(+4.41756058270700, +8.0503895685670, -0.1170591355154);

  // My area Am should coincide with the area below, but my areas Aj, C, and B should be
  // closer to the second configuration than areas below.
  // Am = 17.701, Aj = 26.585, C = 13.819, B = 30.480, 3D: 0.0293438.

  // Compute weights.
  WP2 wachspress;
  std::cout << "2D wp: " << wachspress(query2, vm2, vj2, vp2) << std::endl;
  std::cout << "3D wp: " << wachspress(query3, vm3, vj3, vp3) << std::endl;

  return EXIT_SUCCESS;
}
