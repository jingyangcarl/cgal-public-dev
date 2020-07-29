#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Eigen_diagonalize_traits.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

#if defined(CGAL_USE_OSQP)

// Typedefs.
using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Line_2    = typename Kernel::Line_2;
using Points_2  = std::vector<Point_2>;
using Indices   = std::vector<std::size_t>;
using Segments  = std::vector<Segment_2>;

using NQ = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using AR = CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using OR = CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;

int main(int argc, char *argv[]) {

  // If we want to load a different file, we load it from a path.
  // Each point comes with the index of the corresponding group = region.
  // The file format: x y z i, where i is the group index. The points
  // are 2D hence z = 0. Each group contains points, which form
  // an approximate line.
  std::string path = "data/real_data_2.xyzi";
  if (argc > 1) path = argv[1];
  Saver<Kernel> saver;

  // Initialize input groups = regions.
  std::vector<Points_2> regions;
  initialize_regions(path, regions);

  // Fit a line to each region.
  Line_2 line; Point_2 centroid;
  std::vector<Line_2> lines;
  lines.reserve(regions.size());
  for (const auto& region : regions) {
    CGAL::linear_least_squares_fitting_2(
      region.begin(), region.end(), line, centroid, CGAL::Dimension_tag<0>(),
      Kernel(), CGAL::Eigen_diagonalize_traits<FT, 2>());
    lines.push_back(line);
  }

  // Cut each line at the ends of the corresponding region.
  Segments segments;
  segments.reserve(lines.size());
  Point_2 source, target;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    boundary_points_on_line_2(
      regions[i], lines[i], source, target);
    segments.push_back(Segment_2(source, target));
  }

  // Save input segments.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_real_data_2_before";
    saver.export_eps_segments(segments, full_path, FT(3) / FT(2));
  }

  // Angle regularization.
  const FT max_angle_2 = FT(80);

  // Create neigbor query and angle-based regularization model.
  NQ neighbor_query(segments);
  AR angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2));

  // Regularize.
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);

  std::cout << "* number of modified segments (angles) = " <<
    angle_regularization.number_of_modified_segments() << std::endl;

  // Offset regularization.
  const FT max_offset_2 = FT(2);

  // Get groups of parallel segments after angle regularization.
  std::vector<Indices> pgroups;
  angle_regularization.parallel_groups(
    std::back_inserter(pgroups));

  // Create offset-based regularization model.
  OR offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2));

  // Add each group of parallel segments with at least 2 segments.
  neighbor_query.clear();
  for (const auto& pgroup : pgroups) {
    neighbor_query.add_group(pgroup);
    offset_regularization.add_group(pgroup);
  }

  // Regularize.
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, offset_regularization);

  std::cout << "* number of modified segments (offsets) = " <<
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_real_data_2_after";
    saver.export_eps_segments(segments, full_path, FT(3) / FT(2));
  }
}

#else
int main(void) {
  std::cout << "This example requires the OSQP library." << std::endl;
  return EXIT_SUCCESS;
}
#endif // defined(CGAL_USE_OSQP)
