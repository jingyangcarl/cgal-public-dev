#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Shape_regularization.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Vector_2  = typename Kernel::Vector_2;
using Segment_2 = typename Kernel::Segment_2;
using Line_2    = typename Kernel::Line_2;
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
  CGAL::Shape_regularization::QP_regularization<Kernel, Input_range, Neighbor_query, Angle_regularization, QP_solver>;
using QP_offset_regularizer = 
  CGAL::Shape_regularization::QP_regularization<Kernel, Input_range, Neighbor_query, Offset_regularization, QP_solver>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

void boundary_points_on_line_2(
  const std::vector<Point_2>& points,
  const Line_2& line,
  Point_2 &p, Point_2 &q) {

  const FT max_value = FT(1000000000000);
  
  FT min_proj_value =  max_value;
  FT max_proj_value = -max_value;

  const Vector_2 ref_vector = line.to_vector();
  const Point_2& ref_point  = points[0];
  
  for (const auto& point : points) {
    const Vector_2 curr_vector(ref_point, point);
    const FT value = CGAL::scalar_product(curr_vector, ref_vector);
    if (value < min_proj_value) {
      min_proj_value = value;
      p = point; }
    if (value > max_proj_value) {
      max_proj_value = value;
      q = point; }
  }
}

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize real data 2 example started" 
  << std::endl << std::endl;

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;

  // Load xyz data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/real_data_2.reg");
  CGAL::set_ascii_mode(in);

  if (!in) {
    std::cout << 
    "Error: cannot read the file real_data_2.reg!" << std::endl;
    std::cout << 
    "You can either create a symlink to the data folder or provide this file by hand." 
    << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  Point_2 p; 
  double stub;
  std::size_t region_index;
  std::map<std::size_t, std::vector<Point_2> > point_map;
  while (!in.eof()) {
    in >> p >> stub >> region_index;
    point_map[region_index].push_back(p);
  }
  in.close();

  Line_2 line;
  std::vector<Line_2> lines;
  for (const auto& pair : point_map) {
    const auto& points = pair.second;  
    CGAL::linear_least_squares_fitting_2(
      points.begin(), points.end(), line, CGAL::Dimension_tag<0>());
    lines.push_back(line);
  }

  Point_2 s, t;
  for (std::size_t i = 0; i < lines.size(); ++i) {  
    boundary_points_on_line_2(point_map[i], lines[i], s, t);
    input_range.push_back(Segment_2(s, t));
  }

  // Save input segments.
  if (argc > 2) {
    Saver saver;
    const std::string full_path = 
      std::string(argv[2]) + "regularize_real_data_2_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  // Create a solver.
  QP_solver qp_solver;

  // Create a neighbor query.
  Indices group(input_range.size());
  std::iota(group.begin(), group.end(), 0);
  Neighbor_query neighbor_query(input_range);

  // Angle regularization.
  const FT max_angle = FT(60);
  Angle_regularization angle_regularization(
    input_range, max_angle);

  neighbor_query.add_group(group);
  angle_regularization.add_group(group);

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

  const FT max_distance = FT(95) / FT(100);
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
  if (argc > 2) {
    Saver saver;
    const std::string full_path = 
      std::string(argv[2]) + "regularize_real_data_2_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize real data 2 example finished" 
  << std::endl << std::endl;
}
