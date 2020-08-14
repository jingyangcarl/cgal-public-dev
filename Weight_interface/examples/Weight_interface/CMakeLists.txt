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
    "uniform_weights.cpp")
  create_single_source_cgal_program(
    "three_point_family_weights.cpp")

  create_single_source_cgal_program(
    "inverse_distance_weights.cpp")
  create_single_source_cgal_program(
    "shepard_weights.cpp")

  create_single_source_cgal_program(
    "wachspress_weights.cpp")
  create_single_source_cgal_program(
    "authalic_weights.cpp")

  create_single_source_cgal_program(
    "mean_value_weights.cpp")
  create_single_source_cgal_program(
    "tangent_weights.cpp")

  create_single_source_cgal_program(
    "discrete_harmonic_weights.cpp")
  create_single_source_cgal_program(
    "cotangent_weights.cpp")

  create_single_source_cgal_program(
    "uniform_region_weights.cpp")
  create_single_source_cgal_program(
    "triangular_region_weights.cpp")
  create_single_source_cgal_program(
    "barycentric_region_weights.cpp")
  create_single_source_cgal_program(
    "voronoi_region_weights.cpp")
  create_single_source_cgal_program(
    "mixed_voronoi_region_weights.cpp")

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
