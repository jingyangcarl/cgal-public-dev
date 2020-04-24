#include <vector>
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Segments/Parallel_groups_2.h>

namespace SR = CGAL::Shape_regularization;

template<
typename Saver,
typename Segment_2>
void save_group(
  const std::vector<Segment_2>& segments,
  const std::vector<std::size_t>& group,
  const std::string name,
  Saver& saver) {
  
  std::vector<Segment_2> edges;
  for (const std::size_t seg_index : group)
    edges.push_back(segments[seg_index]);
  saver.export_polylines(edges, "/Users/monet/Documents/gsoc/ggr/logs/" + name);
}

template<class Traits>
void test_parallel_groups() { 
  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  using Segments = std::vector<Segment_2>;
  using PG = SR::Segments::Parallel_groups_2<Traits, Segments>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(1, 1), Point_2(4, 1)), // bottom group
    Segment_2(Point_2(1, 2), Point_2(4, 2)),
    Segment_2(Point_2(1, 3), Point_2(FT(399) / FT(100), FT(319) / FT(100))),

    Segment_2(Point_2(1, 4), Point_2(1, 6)), // top left group
    Segment_2(Point_2(2, 5), Point_2(2, 8)),

    Segment_2(Point_2(3, 5), Point_2(6, 7)), // top right group
    Segment_2(Point_2(7, 6), Point_2(4, 4))
  };
  // saver.export_polylines(segments, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/input");
  
  const PG grouping(
    segments, CGAL::parameters::max_angle(FT(55) / FT(10)));
  std::vector<Indices> groups;
  grouping.parallel_groups(
    std::back_inserter(groups));
  assert(groups.size() == 3);

  // save_group(segments, groups[0], "group0", saver);
  assert(groups[0].size() == 3);
  // save_group(segments, groups[1], "group1", saver);
  assert(groups[1].size() == 2);
  // save_group(segments, groups[2], "group2", saver);
  assert(groups[2].size() == 2);
}

int main() {
  test_parallel_groups< CGAL::Simple_cartesian<double> >();
  test_parallel_groups< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_parallel_groups< CGAL::Exact_predicates_exact_constructions_kernel >();
  return EXIT_SUCCESS;
}
