#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/regularize_segments.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_unique_segments() {

  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Segments  = std::vector<Segment_2>;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  Saver saver;
  const Segments segments = {
    Segment_2(Point_2(1, 1), Point_2(3, 1)), // the right group
    Segment_2(Point_2(4, 1), Point_2(7, 1)),

    Segment_2(Point_2( 2,  2), Point_2( 3,  3)), // the middle group
    Segment_2(Point_2(-1, -1), Point_2(-2, -2)),

    Segment_2(Point_2(-1, 1), Point_2(-1, 3)), // the left group
    Segment_2(Point_2(-FT(3) / FT(2), FT(5) / FT(2)), Point_2(-FT(3) / FT(2), FT(3) / FT(2))),
  };
  // saver.export_eps_segments(segments,
  //   "/Users/monet/Documents/gsoc/ggr/logs/us_input", 100);

  const Segments ref_segments = {
    Segment_2(Point_2( 1,  1), Point_2(7, 1)),
    Segment_2(Point_2(-2, -2), Point_2(3, 3)),
    Segment_2(Point_2(-FT(11) / FT(10), 1), Point_2(-FT(11) / FT(10), 3))
  };

  std::vector<Segment_2> unique;
  SR::Segments::unique_segments(
    segments, std::back_inserter(unique), CGAL::parameters::all_default());
  assert(unique.size() == 4);
  // saver.export_eps_segments(unique, "/Users/monet/Documents/gsoc/ggr/logs/us_output_df", 100);
  assert(unique[0] == ref_segments[0]);
  assert(unique[1] == ref_segments[1]);
  assert(unique[2] == segments[4]);
  assert(unique[3] == segments[5]);

  unique.clear();
  SR::Segments::unique_segments(
    segments, std::back_inserter(unique), CGAL::parameters::max_offset(1));
  assert(unique.size() == 3);
  // saver.export_eps_segments(unique, "/Users/monet/Documents/gsoc/ggr/logs/us_output", 100);
  assert(unique[0] == ref_segments[0]);
  assert(unique[1] == ref_segments[1]);
  assert(unique[2] == ref_segments[2]);

  unique.clear();
  SR::Segments::unique_segments(
    segments, std::back_inserter(unique), CGAL::parameters::
    max_offset(1).
    preserve_order(true));
  assert(unique.size() == 3);
  // saver.export_eps_segments(unique, "/Users/monet/Documents/gsoc/ggr/logs/us_output_pr", 100);
  assert(unique[0] == ref_segments[0]);
  assert(unique[1] == ref_segments[1]);
  assert(unique[2] == ref_segments[2]);
}

int main() {
  test_unique_segments< CGAL::Simple_cartesian<double> >();
  test_unique_segments< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_unique_segments< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_unique_segments: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
