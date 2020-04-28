#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_3_segments() { 

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
    Segment_2(Point_2(0             , 0)              , Point_2(0             , 1)),
    Segment_2(Point_2(FT(1) / FT(10), 0)              , Point_2(FT(2) / FT(10), 1)),
    Segment_2(Point_2(0             , FT(11) / FT(10)), Point_2(FT(2) / FT(10), FT(11) / FT(10)))
  };
  assert(segments.size() == 3);
  // saver.export_polylines(segments, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/3_input");

  NQ neighbor_query(segments, smap);
  neighbor_query.create_unique_group();
  
  const FT max_angle_2 = FT(10);
  AR angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2), smap);
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
  //   "/Users/monet/Documents/gsoc/ggr/logs/3_angles");

  assert(segments.size() == 3);
  assert(parallel_groups.size() == 2);
  assert(orthogonal_groups.size() == 1);
  assert(num_segments_angles == 3);

  std::vector<int> reference_values;
  reference_values.reserve(3);
  reference_values.push_back(1);
  reference_values.push_back(1);
  reference_values.push_back(2);
  SR::Tests::check_reference_values(segments, reference_values);

  const FT max_offset_2 = FT(1) / FT(100);
  OR offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2), smap);

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
  //   "/Users/monet/Documents/gsoc/ggr/logs/3_offsets");

  assert(segments.size() == 3);
  assert(collinear_groups.size() == 3);
  assert(num_segments_offsets == 0);
  SR::Tests::check_reference_values(segments, reference_values);
}

int main() {
  test_3_segments< CGAL::Simple_cartesian<double> >();
  test_3_segments< CGAL::Exact_predicates_inexact_constructions_kernel >(); 
  test_3_segments< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_3_segments: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
