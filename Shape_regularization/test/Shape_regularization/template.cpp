#include "include/Saver.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace SR = CGAL::Shape_regularization;

template<class Traits>
void test_with_name() { 
  using FT        = typename Traits::FT;
  using Point_2   = typename Traits::Point_2;
  using Segment_2 = typename Traits::Segment_2;
  using Indices   = std::vector<std::size_t>;
  using Saver     = SR::Tests::Saver<Traits>;

  // To be added!
}

int main() {
  test_neighbor_query< CGAL::Simple_cartesian<double> >();
  test_neighbor_query< CGAL::Exact_predicates_inexact_constructions_kernel >();
  test_neighbor_query< CGAL::Exact_predicates_exact_constructions_kernel >();

  return EXIT_SUCCESS;
}
