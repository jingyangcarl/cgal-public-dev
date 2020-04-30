#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT        = typename Kernel::FT;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;
using Segments  = std::vector<Segment_2>;

using Neighbor_query = 
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using Angle_regularization = 
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using Offset_regularization = 
  CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;

using Quadratic_program = 
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>;

using QP_angle_regularizer = 
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Angle_regularization, Quadratic_program>;
using QP_offset_regularizer = 
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Offset_regularization, Quadratic_program>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize 15 segments.
  Segments segments;
  CGAL::Shape_regularization::Examples::
  create_example_15(segments);

  // We create three groups of segments: 
  // outer, top and bottom rhombuses.
  std::vector<Indices> groups(3);
  groups[0] = {0, 1, 2, 3, 4, 5, 6}; // outer
  groups[1] = {7, 8, 9, 10};         // top rhombus
  groups[2] = {11, 12, 13, 14};      // bottom rhombus

  // Save input segments.
  if (path != "") {
    const std::string full_path = path + "regularize_15_segments_before";
    saver.export_polylines(segments, full_path);
  }

  // Angle regularization.
  Quadratic_program qp_angles;
  Neighbor_query neighbor_query(segments);

  const FT max_angle_2 = FT(10);
  Angle_regularization angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2));

  for (const auto& group : groups) {
    neighbor_query.add_group(group);
    angle_regularization.add_group(group);
  }

  QP_angle_regularizer qp_angle_regularizer(
    segments, neighbor_query, angle_regularization, qp_angles);
  qp_angle_regularizer.regularize();

  std::cout << "* number of modified segments (angles) = " << 
    angle_regularization.number_of_modified_segments() << std::endl;

  // Offset regularization.
  Quadratic_program qp_offsets;
  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  const FT max_offset_2 = FT(1) / FT(10);
  Offset_regularization offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2));

  neighbor_query.clear();
  for (const auto& group : parallel_groups) {
    neighbor_query.add_group(group);
    offset_regularization.add_group(group);
  }

  QP_offset_regularizer qp_offset_regularizer(
    segments, neighbor_query, offset_regularization, qp_offsets);
  qp_offset_regularizer.regularize();

  std::cout << "* number of modified segments (offsets) = " << 
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (path != "") {
    const std::string full_path = path + "regularize_15_segments_after";
    saver.export_polylines(segments, full_path);
  }
}
