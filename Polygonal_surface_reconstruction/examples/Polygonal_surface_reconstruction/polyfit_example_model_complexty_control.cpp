#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_with_segments.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygonal_surface_reconstruction.h>
#include <CGAL/Timer.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel		Kernel;

typedef CGAL::Point_set_with_segments<Kernel>					Point_set_with_segments;
typedef	CGAL::Polygonal_surface_reconstruction<Kernel>			Polygonal_surface_reconstruction;

typedef Kernel::Point_3											Point;
typedef CGAL::Surface_mesh<Point>								Surface_mesh;


/*
* The following example shows how to control the model complexity by
* increasing the influence of the model complexity term.
* In this example, the intermediate results from plane extraction and 
* candidate generation are cached and reused.
*/

int main()
{
	Point_set_with_segments point_set;

	const std::string& input_file("data/building.vg");
	std::cout << "Loading point cloud: " << input_file << "...";

	CGAL::Timer t;
	t.start();
	if (!point_set.read(input_file)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	} 
	else 
		std::cout << " Done. " 
		<< point_set.number_of_points() << " points, " 
		<< point_set.planar_segments().size() << " planar segments. Time: " << t.time() << " sec." << std::endl;

	Polygonal_surface_reconstruction algo;
	Surface_mesh candidate_faces;
	std::cout << "Generating candidate faces...";
	t.reset();
	if (!algo.generate_candidate_faces(point_set, candidate_faces)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else
		std::cout << " Done. " << candidate_faces.number_of_faces() << " candidate faces. Time: " << t.time() << " sec." << std::endl;

	std::cout << "Computing confidences...";
	t.reset();
	algo.compute_confidences(point_set, candidate_faces);
	std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

	// increasing the weight of the model complexity term results in less detailed 3D models.
	
	// model 1: most details
	Surface_mesh model;
	std::cout << "Optimizing with complexity = 0.2...";
	t.reset();
	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.2)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else {
		const std::string& output_file = "data/building_result_complexity-0.2.off";
		if (CGAL::write_off(std::ofstream(output_file.c_str()), model)) 
			std::cout << " Done. " << model.number_of_faces() << " faces. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl; 
	}

	// model 2: less details
	std::cout << "Optimizing with complexity = 0.4...";
	t.reset();
	if (!algo.select_faces(candidate_faces, model, 0.43, 0.27, 0.4)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else {
		const std::string& output_file = "data/building_result_complexity-0.4.off";
		if (CGAL::write_off(std::ofstream(output_file.c_str()), model))
			std::cout << " Done. " << model.number_of_faces() << " faces. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
	}

	// model 3: least details
	std::cout << "Optimizing with complexity = 0.6...";
	t.reset();
	if (!algo.select_faces(candidate_faces, model, 0.3, 0.1, 0.6)) {
		std::cerr << " Failed." << std::endl;
		return EXIT_FAILURE;
	}
	else {
		const std::string& output_file = "data/building_result_complexity-0.6.off";
		if (CGAL::write_off(std::ofstream(output_file.c_str()), model))
			std::cout << " Done. " << model.number_of_faces() << " faces. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
	}

	return EXIT_SUCCESS;
}
