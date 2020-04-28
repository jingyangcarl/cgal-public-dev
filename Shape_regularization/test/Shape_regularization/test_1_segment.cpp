#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_1_segment() { 

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segments = std::vector<Segment_2>;
  using Segment_map = CGAL::Identity_property_map<Segment_2>;

  using NQ = SR::Segments::Delaunay_neighbor_query_2<Traits, Segments, Segment_map>;
  using AR = SR::Segments::Angle_regularization_2<Traits, Segments, Segment_map>;
  using OR = SR::Segments::Offset_regularization_2<Traits, Segments, Segment_map>;
  using QP = SR::OSQP_quadratic_program<FT>;

  using ARegularizer = SR::QP_regularization<Traits, Segments, NQ, AR, QP>;
  using ORegularizer = SR::QP_regularization<Traits, Segments, NQ, OR, QP>;

  Saver saver;
  Segment_map smap;
  Segments segments = {
    Segment_2(Point_2(1, 1), Point_2(1, 4))
  };
  assert(segments.size() == 1);
  // saver.export_polylines(segments, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/1_input");

  NQ neighbor_query(segments, smap);
  neighbor_query.create_unique_group();
  
  AR angle_regularization(
    segments, CGAL::parameters::all_default(), smap);
  angle_regularization.create_unique_group();

  QP qp_angles;
  ARegularizer aregularizer(
    segments, neighbor_query, angle_regularization, qp_angles);
  aregularizer.regularize();

  std::vector<Indices> parallel_groups;
  angle_regularization.parallel_groups(
    std::back_inserter(parallel_groups));

  std::vector<Indices> orthogonal_groups;
  angle_regularization.orthogonal_groups(
    std::back_inserter(orthogonal_groups));

  const std::size_t num_segments_angles = 
    angle_regularization.number_of_modified_segments();

  // saver.export_polylines(segments, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/1_angles");

  assert(segments.size() == 1);
  assert(parallel_groups.size() == 1);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 0);

  OR offset_regularization(
    segments, CGAL::parameters::all_default(), smap);

  neighbor_query.clear();
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
  //   "/Users/monet/Documents/gsoc/ggr/logs/1_offsets");

  assert(segments.size() == 1);
  assert(collinear_groups.size() == 1);
  assert(num_segments_offsets == 0);
}

int main() {
  test_1_segment< CGAL::Simple_cartesian<double> >();
  test_1_segment< CGAL::Exact_predicates_inexact_constructions_kernel >(); 
  test_1_segment< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_1_segment: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
