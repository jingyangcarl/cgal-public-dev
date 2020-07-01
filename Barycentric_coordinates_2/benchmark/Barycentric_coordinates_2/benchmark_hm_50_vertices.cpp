#include <vector>
#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Delaunay_domain_2.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT       = typename Kernel::FT;
using Point_2  = typename Kernel::Point_2;
using Timer    = CGAL::Real_timer;
using Vertices = std::vector<Point_2>;

using Domain = CGAL::Barycentric_coordinates::Delaunay_domain_2<Vertices, Kernel>;
using HMC2   = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Vertices, Domain, Kernel>;

int main() {

  Timer timer;
  std::cout.precision(10);
  const std::size_t number_of_runs = 1;

  const std::vector<Point_2> vertices = {
    Point_2(0.03, 0.05), Point_2(0.07, 0.04), Point_2(0.10, 0.04),
    Point_2(0.14, 0.04), Point_2(0.17, 0.07), Point_2(0.20, 0.09),
    Point_2(0.22, 0.11), Point_2(0.25, 0.11), Point_2(0.27, 0.10),
    Point_2(0.30, 0.07), Point_2(0.31, 0.04), Point_2(0.34, 0.03),
    Point_2(0.37, 0.02), Point_2(0.40, 0.03), Point_2(0.42, 0.04),
    Point_2(0.44, 0.07), Point_2(0.45, 0.10), Point_2(0.46, 0.13),
    Point_2(0.46, 0.19), Point_2(0.47, 0.26), Point_2(0.47, 0.31),
    Point_2(0.47, 0.35), Point_2(0.45, 0.37), Point_2(0.41, 0.38),
    Point_2(0.38, 0.37), Point_2(0.35, 0.36), Point_2(0.32, 0.35),
    Point_2(0.30, 0.37), Point_2(0.28, 0.39), Point_2(0.25, 0.40),
    Point_2(0.23, 0.39), Point_2(0.21, 0.37), Point_2(0.21, 0.34),
    Point_2(0.23, 0.32), Point_2(0.24, 0.29), Point_2(0.27, 0.24),
    Point_2(0.29, 0.21), Point_2(0.29, 0.18), Point_2(0.26, 0.16),
    Point_2(0.24, 0.17), Point_2(0.23, 0.19), Point_2(0.24, 0.22),
    Point_2(0.24, 0.25), Point_2(0.21, 0.26), Point_2(0.17, 0.26),
    Point_2(0.12, 0.24), Point_2(0.07, 0.20), Point_2(0.03, 0.15),
    Point_2(0.01, 0.10), Point_2(0.02, 0.07)
  };

  std::list<Point_2> list_of_seeds;
  list_of_seeds.push_back(Point_2(0.1, 0.1));

  Domain domain(vertices);
  domain.create(0.001, list_of_seeds);
  std::cout << "benchmark_hm_50_vertices, num vertices: " <<
    domain.number_of_vertices() << std::endl;

  HMC2 harmonic_coordinates_2(vertices, domain);

  timer.start();
  harmonic_coordinates_2.setup();
  timer.stop();
  std::cout << "benchmark_hm_50_vertices, setup (CPU time): " <<
    timer.time() << " seconds" << std::endl;

  timer.reset(); timer.start();
  harmonic_coordinates_2.factorize();
  timer.stop();
  std::cout << "benchmark_hm_50_vertices, factorize (CPU time): " <<
    timer.time() << " seconds" << std::endl;

  timer.reset(); timer.start();
  harmonic_coordinates_2.solve();
  timer.stop();
  std::cout << "benchmark_hm_50_vertices, solve (CPU time): " <<
    timer.time() << " seconds" << std::endl;

  double time = 0.0;
  for (std::size_t i = 0; i < number_of_runs; ++i) {
    harmonic_coordinates_2.clear();
    timer.reset(); timer.start();
    harmonic_coordinates_2.compute();
    timer.stop();
    time += timer.time();
  }

  const double mean_time =
    time / static_cast<double>(number_of_runs);
  std::cout << "benchmark_hm_50_vertices, compute, mean time (CPU time): " <<
    mean_time << " seconds" << std::endl;

  return EXIT_SUCCESS;
}
