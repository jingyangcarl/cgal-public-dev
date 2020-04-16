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

// Choose the type of a solver.
#define USE_OSQP_SOLVER
#if defined(USE_OSQP_SOLVER)
using Quadratic_program = 
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>; // OSQP sparse solver
#else
using Quadratic_program = 
  CGAL::Shape_regularization::CGAL_quadratic_program<FT>; // CGAL dense solver
#endif

using QP_angle_regularizer = 
  CGAL::Shape_regularization::QP_regularization<Kernel, Input_range, Neighbor_query, Angle_regularization, Quadratic_program>;
using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize framework example started" 
  << std::endl << std::endl;

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;

  // Load segments either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/framework.polylines");
  CGAL::set_ascii_mode(in);

  if (!in) {
    std::cout << 
    "Error: cannot read the file framework.polylines!" << std::endl;
    std::cout << 
    "You can either create a symlink to the data folder or provide this file by hand." 
    << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  Point_2 s, t; double stub;
  while (!in.fail()) {
    in >> stub >> s >> stub >> t >> stub;
    input_range.push_back(Segment_2(s, t));
  }
  input_range.erase(input_range.begin() + input_range.size() - 1);
  in.close();

  // Save input segments.
  if (argc > 2) {
    Saver saver;
    const std::string full_path = 
      std::string(argv[2]) + "regularize_framework_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  // Create a solver.
  Quadratic_program qp_angles;

  // Create a neighbor query.
  Indices group(input_range.size());
  std::iota(group.begin(), group.end(), 0);
  Neighbor_query neighbor_query(input_range);

  // Angle regularization.
  const FT max_angle = FT(25);
  Angle_regularization angle_regularization(
    input_range, max_angle);

  neighbor_query.add_group(group);
  angle_regularization.add_group(group);

  QP_angle_regularizer qp_angle_regularizer(
    input_range, neighbor_query, angle_regularization, qp_angles);
  qp_angle_regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of modified segments (angles) = " << 
    angle_regularization.number_of_modified_segments() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized segments.
  if (argc > 2) {
    Saver saver;
    const std::string full_path = 
      std::string(argv[2]) + "regularize_framework_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize framework example finished" 
  << std::endl << std::endl;
}
