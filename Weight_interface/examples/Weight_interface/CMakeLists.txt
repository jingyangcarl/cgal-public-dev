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
    "weight_interface.cpp")
  create_single_source_cgal_program(
    "coordinates_one_query.cpp")
  create_single_source_cgal_program(
    "coordinates_multiple_queries.cpp")
  create_single_source_cgal_program(
    "projection_traits.cpp")
  create_single_source_cgal_program(
    "custom_traits.cpp")
  create_single_source_cgal_program(
    "flattening.cpp")

  find_package(Eigen3 REQUIRED)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})
    create_single_source_cgal_program(
      "weighted_laplacian.cpp")
  else()
    message(WARNING
    "Several examples require the Eigen library, and will not be compiled.")
  endif()

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
