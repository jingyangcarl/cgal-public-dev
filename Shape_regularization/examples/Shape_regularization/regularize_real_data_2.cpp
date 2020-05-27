#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Line_2    = typename Kernel::Line_2;
using Indices   = std::vector<std::size_t>;
using Points_2  = std::vector<Point_2>;

using Segments = std::vector<Segment_2>;
using Segment_map = CGAL::Identity_property_map<Segment_2>;

using Neighbor_query =
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments, Segment_map>;
using Angle_regularization =
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments, Segment_map>;
using Offset_regularization =
  CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments, Segment_map>;

using Quadratic_program =
  CGAL::Shape_regularization::OSQP_quadratic_program<FT>;

using QP_angle_regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Angle_regularization, Quadratic_program>;
using QP_offset_regularizer =
  CGAL::Shape_regularization::QP_regularization<Kernel, Segments, Neighbor_query, Offset_regularization, Quadratic_program>;

using Saver =
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to load a different file, we load it from a path.
  std::string path = "data/real_data_2.rgn";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize input regions.
  Segment_map segment_map;
  std::vector<Points_2> regions;
  CGAL::Shape_regularization::Examples::
  initialize_regions(path, regions);

  // Fit a line to each region.
  Line_2 line;
  std::vector<Line_2> lines;
  lines.reserve(regions.size());
  for (const auto& region : regions) {
    CGAL::linear_least_squares_fitting_2(
      region.begin(), region.end(), line, CGAL::Dimension_tag<0>());
    lines.push_back(line);
  }

  // Cut each line at the ends of the corresponding region.
  Segments segments;
  segments.reserve(lines.size());
  Point_2 source, target;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    CGAL::Shape_regularization::Examples::boundary_points_on_line_2(
      regions[i], lines[i], source, target);
    segments.push_back(Segment_2(source, target));
  }

  // Save input segments.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_real_data_2_before";
    saver.export_eps_segments(segments, full_path, FT(3) / FT(2));
  }

  // Angle regularization.
  Quadratic_program qp_angles;
  Neighbor_query neighbor_query(segments, segment_map);
  neighbor_query.create_unique_group();

  const FT max_angle_2 = FT(80);
  Angle_regularization angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2), segment_map);
  angle_regularization.create_unique_group();

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

  const FT max_offset_2 = FT(2);
  Offset_regularization offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2), segment_map);

  neighbor_query.clear();
  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  QP_offset_regularizer qp_offset_regularizer(
    segments, neighbor_query, offset_regularization, qp_offsets);
  qp_offset_regularizer.regularize();

  std::cout << "* number of modified segments (offsets) = " <<
    offset_regularization.number_of_modified_segments() << std::endl;

  // Save regularized segments.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_real_data_2_after";
    saver.export_eps_segments(segments, full_path, FT(3) / FT(2));
  }
}
