#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization/Contour_regularization_2.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Indices = std::vector<std::size_t>;

using Input_range = std::vector<Point_2>;
namespace SR = CGAL::Shape_regularization;
using Contour_directions_2 = SR::Longest_principal_direction_2<
  Kernel, Input_range>;
using Contour_regularization_2 = SR::Contour_regularization_2<
  Kernel, Input_range, Contour_directions_2, SR::OPEN>;

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
    "regularize open contour example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string in_path = "";
  if (argc > 1) in_path = argv[1];
  const std::string out_path = "/Users/monet/Documents/gsoc/ggr/logs/";

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;
  initialize_segments(in_path, input_range);
  std::cout << "* number of input vertices = " << input_range.size() << std::endl;

  // Save input contour.
  if (out_path != "") {
    Saver saver;
    const std::string full_path = out_path + "regularize_open_contour_before";
    saver.save_open_contour_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  Contour_directions_2 directions(input_range);
  Contour_regularization_2 regularizer(
    input_range, directions, CGAL::parameters::all_default());

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
    const std::string full_path = out_path + "regularize_open_contour_after";
    saver.save_open_contour_2(contour, full_path);
  }

  std::cout << std::endl << 
    "regularize open contour example finished" 
  << std::endl << std::endl;
}
