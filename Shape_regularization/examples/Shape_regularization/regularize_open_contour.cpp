#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Indices = std::vector<std::size_t>;

using Input_range = std::vector<Point_2>;
using Point_map = CGAL::Identity_property_map<Point_2>;

using Contour_directions_2 = 
  CGAL::Shape_regularization::Contours::Longest_direction_2<Kernel, Input_range, Point_map>;
using Contour_regularization_2 = 
  CGAL::Shape_regularization::Contour_regularization_2<Kernel, Input_range, Contour_directions_2, CGAL::Shape_regularization::OPEN, Point_map>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to save the result in a file, we save it in a path.
  std::string in_path = "";
  if (argc > 1) in_path = argv[1];
  const std::string out_path = "/Users/monet/Documents/gsoc/ggr/logs/";

  // Initialize input range.
  Point_map point_map;
  Input_range input_range;
  CGAL::Shape_regularization::Examples::
  initialize_segments<FT, Point_2, Input_range>(in_path, input_range);

  // Save input contour.
  if (out_path != "") {
    Saver saver;
    const std::string full_path = out_path + "regularize_open_contour_before";
    saver.save_open_contour_2(input_range, full_path);
  }

  // Regularize.
  const bool is_closed = false;
  Contour_directions_2 directions(
    input_range, is_closed, point_map);
  Contour_regularization_2 regularizer(
    input_range, directions, CGAL::parameters::all_default(), point_map);

  std::vector<Point_2> contour;
  regularizer.regularize(
    std::back_inserter(contour));

  std::cout << 
    "* number of principal directions = " << 
    directions.number_of_directions() << std::endl;

  // Save regularized contour.
  if (out_path != "") {
    Saver saver;
    const std::string full_path = out_path + "regularize_open_contour_after";
    saver.save_open_contour_2(contour, full_path);
  }
}
