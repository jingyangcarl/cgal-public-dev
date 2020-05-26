#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_open_contour() {

  using FT      = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Saver   = SR::Tests::Saver<Traits>;

  using Contour   = std::vector<Point_2>;
  using Point_map = CGAL::Identity_property_map<Point_2>;

  using CD = SR::Contours::Longest_direction_2<Traits, Contour, Point_map>;
  using OR = SR::Contour_regularization_2<Traits, Contour, CD, SR::OPEN, Point_map>;

  Saver saver;
  Point_map pmap;
  const Contour contour = {
    Point_2(0, 0), Point_2(  FT(5) / FT(10), -FT(1) / FT(20)),
    Point_2(1, 0), Point_2(FT(105) / FT(100), FT(5) / FT(10)),
    Point_2(1, 1), Point_2(0, 1)
  };
  assert(contour.size() == 6);
  // saver.export_open_contour(contour,
  //   "/Users/monet/Documents/gsoc/ggr/logs/op_input");

  const bool is_closed = false;
  CD directions(
    contour, is_closed, pmap);
  OR open_regularizer(
    contour, directions, CGAL::parameters::all_default(), pmap);

  std::vector<Point_2> regularized;
  open_regularizer.regularize(
    std::back_inserter(regularized));
  const std::size_t num_directions =
    directions.number_of_directions();

  // saver.export_open_contour(regularized,
  //   "/Users/monet/Documents/gsoc/ggr/logs/op_output");

  assert(num_directions == 1);
  assert(regularized.size() == 4);
}

int main() {
  test_open_contour< CGAL::Simple_cartesian<double> >();
  test_open_contour< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_open_contour< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_open_contour: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
