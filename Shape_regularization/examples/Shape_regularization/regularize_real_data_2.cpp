#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_regularization.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;

using Input_range = std::vector<Segment_2>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize real data 2 example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;

  /// here

  // Save input segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_real_data_2_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  /// here

  timer.stop();
  std::cout << 
    "* number of modified segments = " << 0 << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_real_data_2_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize real data 2 example finished" 
  << std::endl << std::endl;
}
