#include <vector>
#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

using Kernel  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Contour = std::vector<Point_2>;
using Timer   = CGAL::Real_timer;

using CD = CGAL::Shape_regularization::Contours::Longest_direction_2<Kernel, Contour>;

void benchmark_contours(
  const std::size_t n) {

  FT step = FT(1), d = FT(1) / FT(10);
  FT x = FT(0), y = FT(0);

  Timer timer;
  Contour contour, regularized;
  const FT max_offset_2 = FT(1) / FT(5);

  contour.reserve(n);
  for (std::size_t i = 0; i < n / 2; ++i) {
    x = step * i;
    if (i % 2 == 0) y = FT(0);
    else y = d;
    const Point_2 vertex = Point_2(x, y);
    contour.push_back(vertex);
  }
  contour.push_back(Point_2(x, FT(2)));

  for (std::size_t i = 0; i < n / 2 - 1; ++i) {
    x -= step;
    if (i % 2 == 0) y = FT(2) - d;
    else y = FT(2);
    const Point_2 vertex = Point_2(x, y);
    contour.push_back(vertex);
  }

  regularized.clear();
  timer.start();
  CD closed_directions(contour, true);
  timer.stop();
  const double longest_closed = timer.time();
  timer.reset();

  timer.start();
  CGAL::Shape_regularization::Contours::regularize_closed_contour(
    contour, closed_directions, std::back_inserter(regularized),
    CGAL::parameters::max_offset(max_offset_2));
  timer.stop();
  const double closed_time = timer.time();
  timer.reset();

  regularized.clear();
  timer.start();
  CD open_directions(contour, false);
  timer.stop();
  const double longest_open = timer.time();
  timer.reset();

  timer.start();
  CGAL::Shape_regularization::Contours::regularize_open_contour(
    contour, open_directions, std::back_inserter(regularized),
    CGAL::parameters::max_offset(max_offset_2));
  timer.stop();
  const double open_time = timer.time();
  timer.reset();

  std::cout.precision(10);
  std::cout << "benchmark_contours " << contour.size() << " (CPU time " <<
  "closed/open): " <<
    closed_time << "/" << open_time <<
  " seconds" << std::endl;
}

int main(int argc, char *argv[]) {

  const std::vector<std::size_t> ns = {
    10, 100, 1000, 10000, 100000, 1000000, 10000000
  };
  for (const std::size_t n : ns)
    benchmark_contours(n);
  return EXIT_SUCCESS;
}
