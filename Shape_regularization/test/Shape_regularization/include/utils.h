#ifndef CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H

// STL includes.
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/Random.h>

namespace CGAL {
namespace Shape_regularization {
namespace Tests {

template<typename Segment_2>
void check_reference_values(
  const std::vector<Segment_2>& input_range,
  const std::vector<int>& reference_values) {

  for (std::size_t i = 0; i < input_range.size(); ++i) {
    const auto& segment = input_range[i];

    const auto point1 = segment.source().x() + segment.source().y();
    const auto point2 = segment.target().x() + segment.target().y();

    const int key = static_cast<int>(floor(
      CGAL::to_double(point1 + point2)));
    // std::cout << key << std::endl;
    assert(key == reference_values[i]);
  }
}

double get_coefficient_value(
  const double theta, double& iterator) {

  if (
    theta == 0.0 ||
    theta == CGAL_PI / 2.0 ||
    theta == CGAL_PI ||
    theta == 3.0 * CGAL_PI / 2.0) {

    iterator = 0.0;
  } else if (
    theta == CGAL_PI / 4.0 ||
    theta == 3.0 * CGAL_PI / 4.0 ||
    theta == 5.0 * CGAL_PI / 4.0 ||
    theta == 7.0 * CGAL_PI / 4.0) {

    iterator = 0.22;
  } else if (
    (theta > 0.0 && theta < CGAL_PI / 4.0) ||
    (theta > CGAL_PI / 2.0 && theta < 3.0 * CGAL_PI / 4.0) ||
    (theta > CGAL_PI && theta < 5.0 * CGAL_PI / 4.0) ||
    (theta > 3.0 * CGAL_PI / 2.0 && theta < 7.0 * CGAL_PI / 4.0)) {

    iterator += 0.02;
  } else
    iterator -= 0.02;

  if (theta < CGAL_PI) return -1.0 * iterator;
  return iterator;
}

template<typename Segment_2>
void create_example_offsets(
  std::vector<Segment_2>& segments) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using Point_2 = typename Traits::Point_2;

  segments.clear();
  segments.reserve(100);

  double theta = 0.0, coef = 0.0, iterator = 0.0;
  double theta_step = CGAL_PI / 25.0;

  while (theta < 2.0 * CGAL_PI) {
    const double st = std::sin(theta);
    const double ct = std::cos(theta);

    const Point_2 a = Point_2(0.0, 0.0);
    const Point_2 b = Point_2(ct, st);

    coef = get_coefficient_value(theta, iterator);
    const Point_2 c = Point_2(ct, st + coef);
    const Point_2 d = Point_2(2.0 * ct, 2.0 * st + coef);
    theta += theta_step;

    segments.push_back(Segment_2(a, b));
    segments.push_back(Segment_2(c, d));
  }
}

template<typename Segment_2>
void create_example_angles(
  std::vector<Segment_2>& segments) {

  using Traits = typename Kernel_traits<Segment_2>::Kernel;
  using Point_2 = typename Traits::Point_2;

  using PG = CGAL::Points_on_segment_2<Point_2>;
  using Creator = CGAL::Creator_uniform_2<Point_2, Segment_2>;
  using Segment_iterator = CGAL::Join_input_iterator_2<PG, PG, Creator>;
  using Count_iterator = CGAL::Counting_iterator<Segment_iterator, Segment_2>;

  segments.clear();
  segments.reserve(100);

  // A horizontal like fan.
  PG p1(Point_2(-250,  -50), Point_2(-250,  50), 50);
  PG p2(Point_2( 250, -250), Point_2( 250, 250), 50);

  Segment_iterator t1(p1, p2);
  Count_iterator t1_begin(t1);
  Count_iterator t1_end(t1, 50);
  std::copy(t1_begin, t1_end, std::back_inserter(segments));

  // A vertical like fan.
  PG p3(Point_2( -50, -250), Point_2( 50, -250), 50);
  PG p4(Point_2(-250,  250), Point_2(250,  250), 50);

  Segment_iterator t2(p3, p4);
  Count_iterator t2_begin(t2);
  Count_iterator t2_end(t2, 50);
  std::copy(t2_begin, t2_end, std::back_inserter(segments));
}

} // namespace Tests
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H
