#include <vector>
#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

using Kernel    = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Segments  = std::vector<Segment_2>;
using Indices   = std::vector<std::size_t>;
using Timer     = CGAL::Real_timer;

using NQ = CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Segments>;
using AR = CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Segments>;
using OR = CGAL::Shape_regularization::Segments::Offset_regularization_2<Kernel, Segments>;

void benchmark_qp_segments(
  const std::size_t n,
  const bool regroup) {

  FT step = FT(1), d = FT(1) / FT(10);
  FT x1 = FT(0), y1 = FT(0), x2 = FT(0), y2 = FT(0);
  const std::size_t m = 10;

  Segments segments;
  segments.reserve(n);

  for (std::size_t i = 0; i < n / 2; ++i) {
    x1 = step * i;
    y1 = FT(0);
    if (i % 2 == 0) x2 = step * i;
    else x2 = step * i - d;
    y2 = FT(1);
    const Point_2 source = Point_2(x1, y1);
    const Point_2 target = Point_2(x2, y2);
    segments.push_back(Segment_2(source, target));
  }

  for (std::size_t i = 0; i < n / 2; ++i) {
    x1 = step * i + d;
    y1 = FT(2);
    if (i % 2 == 0) x2 = step * i + d;
    else x2 = step * i;
    y2 = FT(4);
    const Point_2 source = Point_2(x1, y1);
    const Point_2 target = Point_2(x2, y2);
    segments.push_back(Segment_2(source, target));
  }

  Timer timer;
  timer.start();
  NQ neighbor_query(segments);
  if (regroup) {
    Indices group;
    group.reserve(m);
    for (std::size_t i = 0; i < n;) {
      group.clear();
      for (std::size_t j = 0; j < m; ++j)
        group.push_back(i + j);
      neighbor_query.add_group(group);
      i += m;
    }
  }
  timer.stop();
  const double delaunay_time = timer.time();
  timer.reset();

  const FT max_angle_2 = FT(10);
  timer.start();
  AR angle_regularization(
    segments, CGAL::parameters::max_angle(max_angle_2));
  if (regroup) {
    Indices group;
    group.reserve(m);
    for (std::size_t i = 0; i < n;) {
      group.clear();
      for (std::size_t j = 0; j < m; ++j)
        group.push_back(i + j);
      angle_regularization.add_group(group);
      i += m;
    }
  }
  timer.stop();
  const double setup_angle_time = timer.time();
  timer.reset();

  timer.start();
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, angle_regularization);
  timer.stop();
  const double angle_time = timer.time();
  timer.reset();

  const FT max_offset_2 = FT(1) / FT(5);
  timer.start();
  OR offset_regularization(
    segments, CGAL::parameters::max_offset(max_offset_2));
  timer.stop();
  const double setup_offset_time = timer.time();
  timer.reset();

  timer.start();
  std::vector<Indices> pgroups;
  angle_regularization.parallel_groups(
    std::back_inserter(pgroups));
  timer.stop();
  const double init_group_time = timer.time();
  timer.reset();

  timer.start();
  neighbor_query.clear();
  for (const auto& pgroup : pgroups) {
    neighbor_query.add_group(pgroup);
    offset_regularization.add_group(pgroup);
  }
  timer.stop();
  const double add_group_time = timer.time();
  timer.reset();

  timer.start();
  CGAL::Shape_regularization::Segments::regularize_segments(
    segments, neighbor_query, offset_regularization);
  timer.stop();
  const double offset_time = timer.time();
  timer.reset();

  std::cout.precision(10);
  if (regroup)
    std::cout << "grouped: " ;
  std::cout << "benchmark_qp_segments " << segments.size() << " (CPU time " <<
  "delaunay/setup_angles/angles/setup_offsets/offsets): " <<
    delaunay_time << "/" <<
    setup_angle_time << "/" << angle_time << "/" <<
    setup_offset_time << "/" << offset_time <<
  " seconds" << std::endl;
}

int main(int argc, char *argv[]) {

  const std::vector<std::size_t> ns = {
    10, 100, 500, 1000, 5000, 10000, 15000
  };
  std::cout << std::endl;
  for (const std::size_t n : ns) {
    benchmark_qp_segments(n, false);
    benchmark_qp_segments(n, true);
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
