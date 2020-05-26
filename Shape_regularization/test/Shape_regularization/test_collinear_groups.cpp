#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Segments/Collinear_groups_2.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_collinear_groups() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segments = std::vector<Segment_2>;
  using CG = SR::Segments::Collinear_groups_2<Traits, Segments>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(2, -1), Point_2(5, -1)), // bottom group
    Segment_2(Point_2(6, -1), Point_2(7, -1)),

    Segment_2(Point_2(1, 1), Point_2(4, 1)), // middle group
    Segment_2(Point_2(5, FT(11) / FT(10)), Point_2(6, FT(11) / FT(10))),
    Segment_2(Point_2(4, FT(9)  / FT(10)), Point_2(6, FT(9)  / FT(10))),

    Segment_2(Point_2(2, 2), Point_2(3, 3)), // top left group
    Segment_2(Point_2(5, 5), Point_2(4, 4)),

    Segment_2(Point_2(7, 2), Point_2(5, 4)) // top right group
  };
  // saver.export_polylines(segments,
  //   "/Users/monet/Documents/gsoc/ggr/logs/cg_input");

  const CG grouping(
    segments, CGAL::parameters::all_default());
  std::vector<Indices> groups;
  grouping.groups(
    std::back_inserter(groups));
  assert(groups.size() == 4);

  // saver.export_group(segments, groups[0], "/Users/monet/Documents/gsoc/ggr/logs/cg_group0");
  assert(groups[0].size() == 2);
  // saver.export_group(segments, groups[1], "/Users/monet/Documents/gsoc/ggr/logs/cg_group1");
  assert(groups[1].size() == 3);
  // saver.export_group(segments, groups[2], "/Users/monet/Documents/gsoc/ggr/logs/cg_group2");
  assert(groups[2].size() == 2);
  // saver.export_group(segments, groups[3], "/Users/monet/Documents/gsoc/ggr/logs/cg_group3");
  assert(groups[3].size() == 1);
}

int main() {
  test_collinear_groups< CGAL::Simple_cartesian<double> >();
  test_collinear_groups< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_collinear_groups< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_collinear_groups: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
