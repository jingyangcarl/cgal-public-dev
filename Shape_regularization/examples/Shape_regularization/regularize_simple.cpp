#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

using Kernel = CGAL::Simple_cartesian<double>;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;

using NQ = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using AR = CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using QP = CGAL::Shape_regularization::OSQP_quadratic_program<FT>;
using Regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, NQ, AR, QP>;

int main() {

  Segments segments = {
    Segment_2(Point_2(0.0, 0.0), Point_2(1.0, -0.1)),
    Segment_2(Point_2(1.0, 0.0), Point_2(1.0,  1.0)),
    Segment_2(Point_2(1.0, 1.0), Point_2(0.0,  1.1)),
    Segment_2(Point_2(0.0, 1.0), Point_2(0.0,  0.0))
  };

  NQ neighbor_query(segments);
  neighbor_query.create_unique_group();
  AR angle_regularization(
    segments, CGAL::parameters::all_default());
  angle_regularization.create_unique_group();

  QP quadratic_program;
  Regularizer regularizer(
    segments, neighbor_query, angle_regularization, quadratic_program);
  regularizer.regularize();

  std::cout << "* number of modified segments = " <<
    angle_regularization.number_of_modified_segments() << std::endl;
}
