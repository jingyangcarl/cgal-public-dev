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

  const std::vector<Point_2> vertices = {
    Point_2(0, 0), Point_2(1, 0),
    Point_2(1, 1), Point_2(0, 1)
  };

  const std::vector<FT> scales = {
    0.1, 0.08, 0.05, 0.01, 0.008, 0.005, 0.001
  };

  std::list<Point_2> list_of_seeds;
  list_of_seeds.push_back(Point_2(0.5, 0.5));

  Domain domain(vertices);
  HMC2 harmonic_coordinates_2(vertices, domain);

  for (const FT scale : scales) {
    domain.clear();
    domain.create(scale, list_of_seeds);
    std::cout << "benchmark_hm_4_vertices, num vertices: " <<
      domain.number_of_vertices() << std::endl;
    harmonic_coordinates_2.clear();

    timer.reset(); timer.start();
    harmonic_coordinates_2.setup();
    timer.stop();
    const double setup = timer.time();

    timer.reset(); timer.start();
    harmonic_coordinates_2.factorize();
    timer.stop();
    const double factorize = timer.time();

    timer.reset(); timer.start();
    harmonic_coordinates_2.solve();
    timer.stop();
    const double solve = timer.time();

    const double total = setup + factorize + solve;
    std::cout <<
      "benchmark_hm_4_vertices, compute (CPU time setup/factorize/solve/total): " <<
    setup << "/" << factorize << "/" << solve << "/" << total << " seconds" << std::endl;
  }
  return EXIT_SUCCESS;
}
