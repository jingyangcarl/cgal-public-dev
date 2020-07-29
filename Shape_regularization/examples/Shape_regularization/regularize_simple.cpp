#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

#if defined(CGAL_USE_OSQP)

using Kernel    = CGAL::Simple_cartesian<double>;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;

using NQ = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using AR = CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;

int main() {

  Segments segments = {
    Segment_2(Point_2(0.0, 0.0), Point_2(1.0, -0.1)),
    Segment_2(Point_2(1.0, 0.0), Point_2(1.0,  1.0)),
    Segment_2(Point_2(1.0, 1.0), Point_2(0.0,  1.1)),
    Segment_2(Point_2(0.0, 1.0), Point_2(0.0,  0.0))
  };

  NQ neighbor_query(segments);
  AR angle_regularization(
    segments, CGAL::parameters::all_default());
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);

  std::cout << "* number of modified segments = " <<
    angle_regularization.number_of_modified_segments() << std::endl;
}

#else
int main(void) {
  std::cout << "This example requires the OSQP library." << std::endl;
  return EXIT_SUCCESS;
}
#endif // defined(CGAL_USE_OSQP)
