#ifndef CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H
#define CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H

// STL includes.
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/io.h>

namespace CGAL {
namespace Shape_regularization {
namespace Examples {

template<
typename FT,
typename Point_2,
typename Input_range>
void initialize_segments(
  const std::string path,
  Input_range& input_range) {

  input_range.clear();
  std::ifstream in(path.c_str(), std::ios_base::in);
  CGAL::set_ascii_mode(in);
  in.precision(20);

  if (!in) {
    std::cerr << "Error: cannot read the file with data!" << std::endl;
    exit(EXIT_FAILURE);
  }

  double xd, yd;
  std::string tmp;
  while (!in.eof()) {
    in >> tmp >> xd >> yd >> tmp >> tmp >> tmp >> tmp;
    const FT x = static_cast<FT>(xd);
    const FT y = static_cast<FT>(yd);
    const Point_2 point = Point_2(x, y);
    input_range.push_back(point);
  }
  input_range.pop_back();
  in.close();
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

} // namespace Examples
} // namespace Shape_regularization
} // namespace CGAL

#endif // CGAL_SHAPE_REGULARIZATION_EXAMPLES_UTILS_H
