#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

#include <CGAL/Join_input_iterator.h>
#include <CGAL/Counting_iterator.h>
#include <CGAL/function_objects.h>
#include <CGAL/point_generators_2.h>

#include "include/Saver.h"

// Typedefs.
using Kernel = CGAL::Simple_cartesian<double>;

using FT        = typename Kernel::FT;
using Point_2   = typename Kernel::Point_2;
using Segment_2 = typename Kernel::Segment_2;
using Indices   = std::vector<std::size_t>;

using Input_range = std::vector<Segment_2>;

using Neighbor_query = 
  CGAL::Shape_regularization::Segments::Delaunay_neighbor_query_2<Kernel, Input_range>;
using Angle_regularization = 
  CGAL::Shape_regularization::Segments::Angle_regularization_2<Kernel, Input_range>;
using QP_solver = 
  CGAL::Shape_regularization::OSQP_solver<Kernel>;
using QP_angle_regularizer = 
  CGAL::Shape_regularization::QP_regularization
    <Kernel, Input_range, Neighbor_query, Angle_regularization, QP_solver>;

using Saver = 
  CGAL::Shape_regularization::Examples::Saver<Kernel>;

using PG = CGAL::Points_on_segment_2<Point_2>;
using Creator = CGAL::Creator_uniform_2<Point_2, Segment_2>;
using Segment_iterator = CGAL::Join_input_iterator_2<PG, PG, Creator>;
using Count_iterator = CGAL::Counting_iterator<Segment_iterator, Segment_2>;

int main(int argc, char *argv[]) {

  std::cout << std::endl << 
    "regularize 100 segments angles example started" 
  << std::endl << std::endl;

  // If we want to save the result in a file, we save it in a path.
  std::string path = "";
  if (argc > 1) path = argv[1];

  // Initialize a timer.
  CGAL::Timer timer;

  // Initialize input range.
  Input_range input_range;
  input_range.reserve(100);

  // Create a horizontal like fan.
  PG p1(Point_2(-250,  -50), Point_2(-250,  50), 50); // point generator
  PG p2(Point_2( 250, -250), Point_2( 250, 250), 50);
  
  Segment_iterator t1( p1, p2); // segment generator
  
  Count_iterator t1_begin(t1); // count iterator
  Count_iterator t1_end(t1, 50);

  std::copy(t1_begin, t1_end, std::back_inserter(input_range));

  // Create a vertical like fan.
  PG p3(Point_2( -50, -250), Point_2( 50, -250), 50); // point generator
  PG p4(Point_2(-250,  250), Point_2(250,  250), 50);

  Segment_iterator t2( p3, p4); // segment iterator

  Count_iterator t2_begin( t2); // count iterator
  Count_iterator t2_end(t2, 50);

  std::copy(t2_begin, t2_end, std::back_inserter(input_range));

  // Save input segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_100_segments_angles_before";
    saver.save_segments_2(input_range, full_path);
  }

  // Regularize.
  timer.start();

  // Create a solver.
  QP_solver qp_solver;

  // Create a neighbor query.
  Indices group(input_range.size());
  std::iota(group.begin(), group.end(), 0);
  Neighbor_query neighbor_query(input_range);

  // Angle regularization.
  const FT max_angle = FT(40);
  Angle_regularization angle_regularization(
    input_range, max_angle);

  neighbor_query.add_group(group);
  angle_regularization.add_group(group);

  QP_angle_regularizer qp_angle_regularizer(
    input_range, neighbor_query, angle_regularization, qp_solver);
  qp_angle_regularizer.regularize();

  timer.stop();
  std::cout << 
    "* number of modified segments = " << 
    angle_regularization.number_of_modified_segments() << 
    " in time = " << timer.time() << " sec." 
  << std::endl;

  // Save regularized segments.
  if (path != "") {
    Saver saver;
    const std::string full_path = path + "regularize_100_segments_angles_after";
    saver.save_segments_2(input_range, full_path);
  }

  std::cout << std::endl << 
    "regularize 100 segments angles example finished" 
  << std::endl << std::endl;
}
