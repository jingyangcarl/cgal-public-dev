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
    "uniform_weight.cpp")
  create_single_source_cgal_program(
    "wachspress_area_weight.cpp")
  create_single_source_cgal_program(
    "wachspress_cot_weight.cpp")
  create_single_source_cgal_program(
    "authalic_weight.cpp")
  create_single_source_cgal_program(
    "discrete_harmonic_area_weight.cpp")
  create_single_source_cgal_program(
    "discrete_harmonic_cot_weight.cpp")
  create_single_source_cgal_program(
    "cotangent_weight.cpp")
  create_single_source_cgal_program(
    "mean_value_area_weight.cpp")
  create_single_source_cgal_program(
    "mean_value_tan_weight.cpp")
  create_single_source_cgal_program(
    "tangent_weight.cpp")
  create_single_source_cgal_program(
    "area_weight.cpp")

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
