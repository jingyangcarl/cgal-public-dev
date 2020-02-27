#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;

using Input_range = std::vector<Segment_2>;

using Neighbor_query = 
  CGAL::Shape_regularization::Delaunay_neighbor_query_2<Kernel, Input_range>;
using Ordinate_regularization = 
  CGAL::Shape_regularization::Ordinate_regularization_2<Kernel, Input_range>;
using QP_solver = 
  CGAL::Shape_regularization::OSQP_solver<Kernel>;
using QP_ordinate_regularizer = 
  CGAL::Shape_regularization::QP_regularization
    <Kernel, Input_range, Neighbor_query, Ordinate_regularization, QP_solver>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;
using Parallel_groups = 
  CGAL::Shape_regularization::Parallel_groups_2<Kernel, Input_range>;

double get_coef_value(
  const double theta, double& iterator) {
  
  if (
    theta == 0 || 
    theta == CGAL_PI / 2 || 
    theta == CGAL_PI || 
    theta == 3 * CGAL_PI / 2) {
    
    iterator = 0;
  } else if (
    theta == CGAL_PI / 4 || 
    theta == 3 * CGAL_PI / 4 || 
    theta == 5 * CGAL_PI / 4 || 
    theta == 7 * CGAL_PI / 4) {
    
    iterator = 0.22;
  } else if (
    (theta > 0 && theta < CGAL_PI / 4) || 
    (theta > CGAL_PI / 2 && theta < 3 * CGAL_PI / 4) || 
    (theta > CGAL_PI && theta < 5 * CGAL_PI / 4) || 
    (theta > 3 * CGAL_PI / 2 && theta < 7 * CGAL_PI / 4)) {
    
    iterator += 0.02;
  } else
    iterator -= 0.02;

  if (theta < CGAL_PI) return -1 * iterator;
  return iterator;
}

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize 100 segments ordinates example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;
  input_range.reserve(100);

  double theta = 0.0;
  double coef = 0.0;
  double iterator = 0.0;
  double theta_step = CGAL_PI / 25.0;

  while (theta < 2 * CGAL_PI) {
    const double st = std::sin(theta);
    const double ct = std::cos(theta);
    const Point_2 a = Point_2(0.0, 0.0);
    const Point_2 b = Point_2(ct, st);

    coef = get_coef_value(theta, iterator);
    const Point_2 c = Point_2(ct, st + coef);
    const Point_2 d = Point_2(2 * ct, 2 * st + coef);
    theta += theta_step;

    input_range.push_back(Segment_2(a, b));
    input_range.push_back(Segment_2(c, d));
  }

  // Save input segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_100_segments_ordinates_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  // Create parallel groups.
  const FT tolerance = FT(1);
  Parallel_groups grouping(
    input_range, tolerance);

  std::vector<Indices> parallel_groups;
  grouping.parallel_groups(
    std::back_inserter(parallel_groups));

  std::cout << 
    "* number of parallel groups = " << parallel_groups.size() 
  << std::endl;

  // Create a solver.
  QP_solver qp_solver;

  // Create a neighbor query.
  Neighbor_query neighbor_query(input_range);

  // Ordinate regularization.
  const FT max_distance = FT(25) / FT(100);
  Ordinate_regularization ordinate_regularization(
    input_range, max_distance);

  for (const auto& group : parallel_groups) {
    if (group.size() < 2) continue;
    neighbor_query.add_group(group);
    ordinate_regularization.add_group(group);
  }

  QP_ordinate_regularizer qp_ordinate_regularizer(
    input_range, neighbor_query, ordinate_regularization, qp_solver);
  qp_ordinate_regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of modified segments = " << 
    ordinate_regularization.number_of_modified_segments() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_100_segments_ordinates_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize 100 segments ordinates example finished" 
  << std::endl << std::endl;
}
