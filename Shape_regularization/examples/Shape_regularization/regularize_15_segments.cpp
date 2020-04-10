#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;

using Input_range = std::vector<Segment_2>;

using Neighbor_query = 
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Input_range>;
using Angle_regularization = 
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Input_range>;
using Offset_regularization = 
  CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Input_range>;

using QP_solver = 
  CGAL::Shape_regularization::OSQP_solver<Kernel>;

using QP_angle_regularizer = 
  CGAL::Shape_regularization::QP_regularization
    <Kernel, Input_range, Neighbor_query, Angle_regularization, QP_solver>;
using QP_offset_regularizer = 
  CGAL::Shape_regularization::QP_regularization
    <Kernel, Input_range, Neighbor_query, Offset_regularization, QP_solver>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize 15 segments example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.

  // Initialize 30 points.
  std::vector<Point_2> points = {
    Point_2(1.000000, 1.000000), Point_2(0.925377, 2.995179),
    Point_2(1.000000, 3.000000), Point_2(1.066662, 4.951894),
    Point_2(1.000000, 5.000000), Point_2(2.950000, 4.930389),
    Point_2(3.000000, 4.950000), Point_2(2.934996, 3.008203),
    Point_2(3.085452, 3.003266), Point_2(2.969782, 1.002004),
    Point_2(0.948866, 3.033161), Point_2(2.900000, 3.000000),
    Point_2(0.930000, 1.000000), Point_2(2.860000, 1.002004),
    Point_2(1.600000, 4.000000), Point_2(1.932136, 4.364718),
    Point_2(1.598613, 3.982686), Point_2(2.018220, 3.686595),
    Point_2(1.951872, 4.363094), Point_2(2.290848, 4.054154),
    Point_2(2.018220, 3.686595), Point_2(2.304517, 4.045054),
    Point_2(1.642059, 1.928505), Point_2(1.993860, 2.247986),
    Point_2(1.993860, 2.247986), Point_2(2.259099, 1.919966),
    Point_2(1.629845, 1.923077), Point_2(1.968759, 1.599174),
    Point_2(2.259099, 1.919966), Point_2(1.968759, 1.599174)
  };

  // Create 15 segments.
  Input_range input_range = {
    Segment_2(points[0] , points[1]) , Segment_2(points[2] , points[3]) ,
    Segment_2(points[4] , points[5]) , Segment_2(points[6] , points[7]) ,
    Segment_2(points[8] , points[9]) , Segment_2(points[10], points[11]),
    Segment_2(points[12], points[13]), Segment_2(points[14], points[15]),
    Segment_2(points[16], points[17]), Segment_2(points[18], points[19]),
    Segment_2(points[20], points[21]), Segment_2(points[22], points[23]),
    Segment_2(points[24], points[25]), Segment_2(points[26], points[27]),
    Segment_2(points[28], points[29])
  };

  // We also create three groups of segments: outer, top and bottom rhombuses.
  std::vector<Indices> groups(3);
  groups[0] = {0, 1, 2, 3, 4, 5, 6}; // outer
  groups[1] = {7, 8, 9, 10}; // top rhombus
  groups[2] = {11, 12, 13, 14}; // bottom rhombus

  // Save input segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_15_segments_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  // Create a solver.
  QP_solver qp_solver;

  // Create a neighbor query.
  Neighbor_query neighbor_query(input_range);

  // Angle regularization.
  const FT max_angle = FT(385) / FT(100);
  Angle_regularization angle_regularization(
    input_range, max_angle);

  for (const auto& group : groups) {
    neighbor_query.add_group(group);
    angle_regularization.add_group(group);
  }

  QP_angle_regularizer qp_angle_regularizer(
    input_range, neighbor_query, angle_regularization, qp_solver);
  qp_angle_regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of modified segments (angles) = " << 
    angle_regularization.number_of_modified_segments() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Offset regularization.
  timer.reset(); timer.start();
  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  std::cout << 
    "* number of parallel_groups = " << parallel_groups.size() 
  << std::endl;

  const FT max_distance = FT(1) / FT(10);
  Offset_regularization offset_regularization(
    input_range, max_distance);

  neighbor_query.clear();
  for (const auto& group : parallel_groups) {
    if (group.size() < 2) continue;
    neighbor_query.add_group(group);
    offset_regularization.add_group(group);
  }

  QP_offset_regularizer qp_offset_regularizer(
    input_range, neighbor_query, offset_regularization, qp_solver);
  qp_offset_regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of modified segments (offsets) = " << 
    offset_regularization.number_of_modified_segments() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_15_segments_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize 15 segments example finished" 
  << std::endl << std::endl;
}
