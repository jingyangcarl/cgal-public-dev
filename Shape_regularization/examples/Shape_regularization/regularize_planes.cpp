#include "include/Saver.h"
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;
using Plane_3 = typename Kernel::Plane_3;

using Point_map = CGAL::Identity_property_map<Point_3>;
using Plane_map = CGAL::Identity_property_map<Plane_3>;

using Point_range = std::vector<Point_3>;
using Plane_range = std::vector<Plane_3>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  // Initialize input range.
  Point_range points;
  Plane_range planes;

  // To be added later!

  // Regularize.
  CGAL::Shape_regularization::Planes::regularize_planes(
    points, planes,
    CGAL::parameters::
    point_map(Point_map()).
    plane_map(Plane_map()).
    // point_to_plane_index_map(
    //   CGAL::Shape_regularization::Planes::Point_to_shape_index_map<Kernel>(
    //     points, planes)). // Do we need it at all?
    regularize_parallelism(true).
    regularize_orthogonality(true).
    regularize_coplanarity(false).
    regularize_z_symmetry(true).
    max_angle(FT(10)));
}
