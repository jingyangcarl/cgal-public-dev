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

using Parallel_groups =
  CGAL::Shape_regularization::Segments::Parallel_groups_2<Kernel, Segments>;

using Neighbor_query =
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using Offset_regularization =
  CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;
using Quadratic_program =
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>;
using QP_offset_regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Offset_regularization, Quadratic_program>;

using Saver =
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize 100 segments in a circle.
  Segments segments;
  CGAL::Shape_regularization::Examples::
  create_example_offsets(segments);

  // Save input segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_offsets_before";
    saver.export_eps_segments(segments, full_path);
  }

  // Create parallel groups.
  const FT max_angle_2 = FT(1);
  const Parallel_groups grouping(segments,
  CGAL::parameters::max_angle(max_angle_2));
  std::vector<Indices> parallel_groups;
  grouping.groups(
    std::back_inserter(parallel_groups));

  // Create a solver.
  Quadratic_program qp_offsets;

  // Create a neighbor query.
  Neighbor_query neighbor_query(segments);

  // Offset regularization.
  const FT max_offset_2 = FT(1) / FT(4);
  Offset_regularization offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2));

  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  QP_offset_regularizer qp_offset_regularizer(
    segments, neighbor_query, offset_regularization, qp_offsets);
  qp_offset_regularizer.regularize();

  std::cout << "* number of modified segments = " <<
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_offsets_after";
    saver.export_eps_segments(segments, full_path);
  }
}
