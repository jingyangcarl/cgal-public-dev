# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Barycentric_coordinates_2_Tests)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program(
    "test_almost_degenerate_segment.cpp")
  create_single_source_cgal_program(
    "test_almost_degenerate_triangle.cpp")
  create_single_source_cgal_program(
    "test_boundary_coordinates_at_vertices.cpp")
  create_single_source_cgal_program(
    "test_boundary_coordinates_on_edges.cpp")
  create_single_source_cgal_program(
    "test_dh_almost_degenerate_polygon.cpp")
  create_single_source_cgal_program(
    "test_dh_const_linear_precision.cpp")
  create_single_source_cgal_program(
    "test_dh_triangle.cpp")
  create_single_source_cgal_program(
    "test_dh_weights.cpp")
  create_single_source_cgal_program(
    "test_mv_const_linear_precision.cpp")
  create_single_source_cgal_program(
    "test_mv_special_points.cpp")
  create_single_source_cgal_program(
    "test_mv_triangle.cpp")
  create_single_source_cgal_program(
    "test_mv_weakly_convex_polygon.cpp")
  create_single_source_cgal_program(
    "test_mv_weights.cpp")
  create_single_source_cgal_program(
    "test_segment_coordinates.cpp")
  create_single_source_cgal_program(
    "test_segment_coordinates_with_offset.cpp")
  create_single_source_cgal_program(
    "test_triangle_coordinates.cpp")
  create_single_source_cgal_program(
    "test_wp_almost_degenerate_polygon.cpp")
  create_single_source_cgal_program(
    "test_wp_const_linear_precision.cpp")
  create_single_source_cgal_program(
    "test_wp_dh_unit_square.cpp")
  create_single_source_cgal_program(
    "test_wp_triangle.cpp")
  create_single_source_cgal_program(
    "test_wp_weights.cpp")

else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
