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

using Contour_regularization_2 = CGAL::Shape_regularization::
  Contour_regularization_2<Kernel, Input_range, Point_map>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize closed contour example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;
  Point_map point_map;

  

  // Save input contour.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_closed_contour_before";
    saver.save_closed_contour_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  Contour_regularization_2 regularizer(
    input_range, point_map, true);
  regularizer.estimate_principal_directions();
  regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of principal directions = " << 
    regularizer.number_of_principal_directions() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized contour.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_closed_contour_after";
    saver.save_closed_contour_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize closed contour example finished" 
  << std::endl << std::endl;
}
