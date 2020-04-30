#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT        = typename Kernel::FT;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;
using Segments  = std::vector<Segment_2>;

using Neighbor_query = 
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using Angle_regularization = 
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using Quadratic_program = 
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>;
using QP_angle_regularizer = 
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Angle_regularization, Quadratic_program>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize 100 near-orthogonal segments.
  Segments segments;
  CGAL::Shape_regularization::Examples::
  create_example_angles(segments);

  // Save input segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_angles_before";
    saver.export_polylines(segments, full_path);
  }

  // Create a solver.
  Quadratic_program qp_angles;

  // Create a neighbor query.
  Neighbor_query neighbor_query(segments);
  neighbor_query.create_unique_group();

  // Angle regularization.
  const FT max_angle_2 = FT(40);
  Angle_regularization angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2));
  angle_regularization.create_unique_group();

  QP_angle_regularizer qp_angle_regularizer(
    segments, neighbor_query, angle_regularization, qp_angles);
  qp_angle_regularizer.regularize();

  std::cout << "* number of modified segments = " << 
    angle_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_angles_after";
    saver.export_polylines(segments, full_path);
  }
}
