#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

using Kernel    = CGAL::Simple_cartesian<double>;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;

int main() {

  // Create input segments.
  std::vector<Segment_2> segments = {
    Segment_2(Point_2(0.0, 0.0), Point_2(1.0, -0.1)),
    Segment_2(Point_2(1.0, 0.0), Point_2(1.0,  1.0)),
    Segment_2(Point_2(1.0, 1.0), Point_2(0.0,  1.1)),
    Segment_2(Point_2(0.0, 1.0), Point_2(0.0,  0.0))
  };
  auto copied = segments;
  std::vector<Segment_2> parallel = {
    Segment_2(Point_2(0.0, 0.0), Point_2(1.0, 0.0)),
    Segment_2(Point_2(0.0, 0.1), Point_2(1.0, 0.1))
  };

  // Regularize all segments: both angles and offsets.
  CGAL::Shape_regularization::Segments::regularize_segments(segments);

  // Regularize only angles.
  CGAL::Shape_regularization::Segments::regularize_angles(copied);

  // Regularize only offsets. All input segments must be parallel!
  CGAL::Shape_regularization::Segments::regularize_offsets(parallel);
}
