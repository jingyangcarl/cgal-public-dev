# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Weight_interface_Examples)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "uniform_weight_2.cpp")
  create_single_source_cgal_program(
    "three_point_family_weight_2.cpp")

  create_single_source_cgal_program(
    "inverse_distance_weight_2.cpp")
  create_single_source_cgal_program(
    "shepard_weight_2.cpp")

  create_single_source_cgal_program(
    "wachspress_weight_2.cpp")
  create_single_source_cgal_program(
    "authalic_weight_2.cpp")

  create_single_source_cgal_program(
    "mean_value_weight_2.cpp")
  create_single_source_cgal_program(
    "tangent_weight_2.cpp")

  create_single_source_cgal_program(
    "discrete_harmonic_weight_2.cpp")
  create_single_source_cgal_program(
    "cotangent_weight_2.cpp")

  create_single_source_cgal_program(
    "uniform_cell_weight_2.cpp")
  create_single_source_cgal_program(
    "triangle_cell_weight_2.cpp")
  create_single_source_cgal_program(
    "barycentric_cell_weight_2.cpp")
  create_single_source_cgal_program(
    "voronoi_cell_weight_2.cpp")
  create_single_source_cgal_program(
    "mixed_voronoi_cell_weight_2.cpp")

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
