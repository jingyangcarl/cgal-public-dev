#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Contours/User_defined_directions_2.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_user_defined_directions() { 

  using FT      = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Direction_2 = typename Traits::Direction_2;
  using Saver   = SR::Tests::Saver<Traits>;

  using Contour = std::vector<Point_2>;
  using Point_map = CGAL::Identity_property_map<Point_2>;
  using UD = SR::Contours::User_defined_directions_2<Traits, Contour, Point_map>;
  
  Saver saver;
  Point_map pmap;
  const Contour contour = {
    Point_2( 1, 1), Point_2(4, 1), 
    Point_2( 4, 4), Point_2(7, 1),
    Point_2(10, 4), Point_2(7, 7),
    Point_2(1, 7)
  };
  assert(contour.size() == 7);
  // saver.export_closed_contour(contour, 
  //   "/Users/monet/Documents/gsoc/ggr/logs/ud_input");

  const std::vector<Direction_2> dirs = {
    Direction_2(1, 0),
    Direction_2(-FT(7) / FT(10), FT(7) / FT(10))
  };

  const bool is_closed = true;
  UD closed_directions(
    dirs, contour,  is_closed, pmap);
  UD open_directions(
    dirs, contour, !is_closed, pmap);

  const std::size_t num_closed_directions = 
    closed_directions.number_of_directions();
  const std::size_t num_open_directions = 
    open_directions.number_of_directions();

  assert(num_closed_directions == 2);
  assert(num_closed_directions == num_open_directions);
  assert(num_closed_directions == dirs.size());

  const auto& closed_dirs = closed_directions.get_directions();
  const auto& open_dirs = open_directions.get_directions();

  assert(closed_dirs[0] == open_dirs[0]);
  assert(closed_dirs[1] == open_dirs[1]);
}

int main() {
  test_user_defined_directions< CGAL::Simple_cartesian<double> >();
  test_user_defined_directions< CGAL::Exact_predicates_inexact_constructions_kernel >(); 
  test_user_defined_directions< CGAL::Exact_predicates_exact_constructions_kernel >();
  std::cout << "test_user_defined_directions: SUCCESS" << std::endl;
  return EXIT_SUCCESS;
}
