#include "include/utils.h"
#include "include/Saver.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_2 = typename Kernel::Point_2;
using Contour = std::vector<Point_2>;

using Point_map = CGAL::Identity_property_map<Point_2>;
using Contour_directions = 
  CGAL::Shape_regularization::Contours::Longest_direction_2<Kernel, Contour, Point_map>;
using Contour_regularizer = 
  CGAL::Shape_regularization::Contour_regularization_2<Kernel, Contour, Contour_directions, 
  CGAL::Shape_regularization::OPEN, Point_map>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // If we want to load a different file, we load it from a path.
  std::string path = "data/contour.polylines";
  if (argc > 1) path = argv[1];
  Saver saver;

  // Initialize contour.
  Contour contour;
  CGAL::Shape_regularization::Examples::
  initialize_contour(path, contour);

  // Save input contour.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_open_contour_before";
    saver.export_open_contour(contour, full_path);
  }

  // Regularize.
  const bool is_closed = false;
  Contour_directions directions(
    contour, is_closed, Point_map());

  Contour_regularizer regularizer(
    contour, directions, CGAL::parameters::all_default(), Point_map());

  Contour regularized;
  regularizer.regularize(
    std::back_inserter(regularized));

  std::cout << "* number of directions = " << 
    directions.number_of_directions() << std::endl;

  // Save regularized contour.
  if (argc > 2) {
    const std::string full_path = std::string(argv[2]) + "regularize_open_contour_after";
    saver.export_open_contour(regularized, full_path);
  }
}
