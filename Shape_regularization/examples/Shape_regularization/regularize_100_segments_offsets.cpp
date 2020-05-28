#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

// Typedefs.
using Kernel    = CGAL::Simple_cartesian<double>;
using FT        = typename Kernel::FT;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;
using Indices   = std::vector<std::size_t>;

using PG    = CGAL::Shape_regularization::Segments::Parallel_groups_2<Kernel, Segments>;
using NQ    = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using OR    = CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;
using Saver = CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize 100 segments in a fuzzy circle.
  Segments segments;
  CGAL::Shape_regularization::Examples::
  create_example_offsets(segments);

  // Save input segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_offsets_before";
    saver.export_eps_segments(segments, full_path, FT(100));
  }

  // Find groups of parallel segments.
  const FT max_angle_2 = FT(1);

  const PG grouping(segments,
  CGAL::parameters::max_angle(max_angle_2));
  std::vector<Indices> pgroups;
  grouping.groups(
    std::back_inserter(pgroups));

  // Offset regularization.
  const FT max_offset_2 = FT(1) / FT(4);

  // Create neigbor query and offset-based regularization model.
  NQ neighbor_query(segments);
  OR offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2));

  // Add each group of parallel segments with at least 2 segments.
  for (const auto& pgroup : pgroups) {
    neighbor_query.add_group(pgroup);
    offset_regularization.add_group(pgroup);
  }

  // Regularize.
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, offset_regularization);

  std::cout << "* number of modified segments = " <<
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (path != "") {
    const std::string full_path = path + "regularize_100_segments_offsets_after";
    saver.export_eps_segments(segments, full_path, FT(100));
  }
}
