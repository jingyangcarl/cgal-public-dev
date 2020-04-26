#include <vector>
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Segments/Orthogonal_groups_2.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_orthogonal_groups() { 
  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segments = std::vector<Segment_2>;
  using OG = SR::Segments::Orthogonal_groups_2<Traits, Segments>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(2, 4), Point_2(2, 1)), // left group
    Segment_2(Point_2(4, 1), Point_2(FT(39) / FT(10), 4)),
    Segment_2(Point_2(2, 6), Point_2(7, 6)),
    Segment_2(Point_2(5, 5), Point_2(5, 4)),
    Segment_2(Point_2(5, 3), Point_2(7, 3)),

    Segment_2(Point_2(7, 1), Point_2(10, 0)), // bottom right group

    Segment_2(Point_2( 8, 3), Point_2(FT(98) / FT(10), FT(52) / FT(10))), // top right group
    Segment_2(Point_2(12, 2), Point_2(10, 4))
  };
  // saver.export_polylines(segments, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/input");
  
  const OG grouping(
    segments, CGAL::parameters::all_default());
  std::vector<Indices> groups;
  grouping.orthogonal_groups(
    std::back_inserter(groups));
  assert(groups.size() == 3);

  // saver.export_group(segments, groups[0], "group0");
  assert(groups[0].size() == 5);
  // saver.export_group(segments, groups[1], "group1");
  assert(groups[1].size() == 2);
  // saver.export_group(segments, groups[2], "group2");
  assert(groups[2].size() == 1);
}

int main() {
  test_orthogonal_groups< CGAL::Simple_cartesian<double> >();
  test_orthogonal_groups< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_orthogonal_groups< CGAL::Exact_predicates_exact_constructions_kernel >();
  return EXIT_SUCCESS;
}
