#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

using Kernel    = CGAL::Simple_cartesian<double>;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;

using Neighbor_query =
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using Angle_regularization =
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;

int main() {

  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0.0, 0.0), Point_2(1.0, -0.1)),
    Segment_2(Point_2(1.0, 0.0), Point_2(1.0,  1.0)),
    Segment_2(Point_2(1.0, 1.0), Point_2(0.0,  1.1)),
    Segment_2(Point_2(0.0, 1.0), Point_2(0.0,  0.0))
  };

  Neighbor_query neighbor_query(segments);
  Angle_regularization angle_regularization(
    segments, CGAL::parameters::all_default());
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);

  std::cout << "* number of modified segments = " <<
    angle_regularization.number_of_modified_segments() << std::endl;
}
