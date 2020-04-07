#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Indices = std::vector<std::size_t>;

using Input_range = std::vector<Point_2>;
using Point_map = CGAL::Identity_property_map<Point_2>;

namespace SR = CGAL::Shape_regularization;
using Contour_regularization_2 = SR::Contour_regularization_2<
  Kernel, Input_range, SR::CLOSED>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

// TODO:
// 1. Clean up this example, e.g. remove out_path.

void initialize_segments(
  const std::string path,
  Input_range& input_range) {

  input_range.clear();
  std::ifstream in(path.c_str(), std::ios_base::in);
  CGAL::set_ascii_mode(in);
  in.precision(20);

  if (!in) {
    std::cerr << "Error: Error loading file with data!" << std::endl;
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

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize closed contour example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string in_path = "";
  if (argc > 1) in_path = argv[1];
  const std::string out_path = "/Users/monet/Documents/gsoc/ggr/logs/";
  
  // Parameters.
  const FT min_length_2 = FT(2);
  const FT max_angle_2  = FT(20);
  const FT max_offset_2 = FT(1) / FT(2);

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;
  Point_map point_map;
  initialize_segments(in_path, input_range);
  std::cout << "* number of input vertices = " << input_range.size() << std::endl;

  // Save input contour.
  if (out_path != "") {
    Saver saver;
    const std::string full_path = out_path + "regularize_closed_contour_before";
    saver.save_closed_contour_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  Contour_regularization_2 regularizer(
    input_range, 
  CGAL::parameters::
  point_map(point_map).
  min_length(min_length_2).
  max_angle(max_angle_2).
  max_offset(max_offset_2));

  regularizer.estimate_principal_directions(
    SR::Direction_type::LONGEST);

  std::vector<Point_2> contour;
  regularizer.regularize(
    std::back_inserter(contour));

  timer.stop();
  std::cout << 
    "* number of principal directions = " << 
    regularizer.number_of_principal_directions() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized contour.
  if (out_path != "") {
    Saver saver;
    const std::string full_path = out_path + "regularize_closed_contour_after";
    saver.save_closed_contour_2(contour, full_path);
  }

  std::cout << std::endl << 
    "regularize closed contour example finished" 
  << std::endl << std::endl;
}
