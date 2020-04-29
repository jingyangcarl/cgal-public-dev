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
#include <CGAL/Aff_transformation_3.h>
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

// Change this function!
double get_coef_value(
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

} // namespace Tests
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_TESTS_UTILS_H
