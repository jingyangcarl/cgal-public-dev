#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_100_segments_offsets() {

  using FT        = typename Traits::FT;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segments = std::vector<Segment_2>;
  using Segment_map = CGAL::Identity_property_map<Segment_2>;

  using PG = SR::Segments::Parallel_groups_2<Traits, Segments, Segment_map>;
  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments, Segment_map>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments, Segment_map>;
  using QP = SR::OSQP_quadratic_program<FT>;

  using ORegularizer = SR::QP_regularization<Traits, Segments, NQ, OR, QP>;

  Saver saver;
  Segment_map smap;
  Segments segments;
  SR::Tests::create_example_offsets(segments);
  assert(segments.size() == 100);
  // saver.export_polylines(segments,
  //   "/Users/monet/Documents/gsoc/ggr/logs/100o_input");

  const FT max_angle_2 = FT(1);
  PG grouping(
    segments,
    CGAL::parameters::max_angle(max_angle_2),
    smap);

  std::vector<Indices> parallel_groups;
  grouping.groups(std::back_inserter(
    parallel_groups));

  // Segments output;
  // for (std::size_t i = 0; i < parallel_groups.size(); ++i) {
  //   const auto& parallel_group = parallel_groups[i];
  //   output.clear();
  //   for (const std::size_t idx : parallel_group)
  //     output.push_back(segments[idx]);
  //   saver.export_polylines(output,
  //   "/Users/monet/Documents/gsoc/ggr/logs/output_" + std::to_string(i));
  // }

  assert(segments.size() == 100);
  assert(parallel_groups.size() == 25);

  const FT max_offset_2 = FT(1) / FT(4);
  OR offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2), smap);

  NQ neighbor_query(segments, smap);
  for (const auto& parallel_group : parallel_groups) {
    neighbor_query.add_group(parallel_group);
    offset_regularization.add_group(parallel_group);
  }

  QP qp_offsets;
  ORegularizer oregularizer(
    segments, neighbor_query, offset_regularization, qp_offsets);
  oregularizer.regularize();

  std::vector<Indices> collinear_groups;
  offset_regularization.collinear_groups(
    std::back_inserter(collinear_groups));

  const std::size_t num_segments_offsets =
    offset_regularization.number_of_modified_segments();

  // saver.export_polylines(segments,
  //   "/Users/monet/Documents/gsoc/ggr/logs/100o_offsets");

  assert(segments.size() == 100);
  assert(collinear_groups.size() == 25); // should be the same as parallel_groups!
  assert(num_segments_offsets == 100);
}

int main() {
  test_100_segments_offsets< CGAL::Simple_cartesian<double> >();
  test_100_segments_offsets< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_100_segments_offsets< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_100_segments_offsets: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
