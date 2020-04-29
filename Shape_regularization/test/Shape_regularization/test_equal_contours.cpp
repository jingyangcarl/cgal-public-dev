#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_equal_contours() { 

  using FT      = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Saver   = SR::Tests::Saver<Traits>;

  using Contour   = std::vector<Point_2>;
  using Point_map = CGAL::Identity_property_map<Point_2>;

  using CD = SR::Contours::Longest_direction_2<Traits, Contour, Point_map>;
  using OD = SR::Contours::Longest_direction_2<Traits, Contour, Point_map>;

  using CR = SR::Contour_regularization_2<Traits, Contour, CD, SR::CLOSED, Point_map>;
  using OR = SR::Contour_regularization_2<Traits, Contour, OD, SR::OPEN, Point_map>;
  
  Saver saver;
  Point_map pmap;
  const Contour contour = {
    Point_2(0, 0), Point_2(  FT(5) / FT(10), -FT(1) / FT(20)),
    Point_2(1, 0), Point_2(FT(105) / FT(100), FT(5) / FT(10)),
    Point_2(1, 1), Point_2(0, 1)
  };
  assert(contour.size() == 6);

  const bool is_closed = true;
  CD closed_directions(
    contour,  is_closed, pmap);
  OD open_directions(
    contour, !is_closed, pmap);
    
  CR closed_regularizer(
    contour, closed_directions, CGAL::parameters::all_default(), pmap);
  OR open_regularizer(
    contour, open_directions, CGAL::parameters::all_default(), pmap);

  std::vector<Point_2> closed_contour;
  closed_regularizer.regularize(
    std::back_inserter(closed_contour));

  std::vector<Point_2> open_contour;
  open_regularizer.regularize(
    std::back_inserter(open_contour));
  
  assert(closed_contour.size() == 4);
  assert(closed_contour.size() == open_contour.size());
  for (std::size_t i = 1; i < closed_contour.size(); ++i) {
    const auto& p = closed_contour[i];
    const auto& q = open_contour[i];

    const FT sq_dist = CGAL::squared_distance(p, q);
    assert(sq_dist <= FT(1) / FT(100000));
  }
}

int main() {
  test_equal_contours< CGAL::Simple_cartesian<double> >();
  test_equal_contours< CGAL::Exact_predicates_inexact_constructions_kernel >(); 
  test_equal_contours< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_equal_contours: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
