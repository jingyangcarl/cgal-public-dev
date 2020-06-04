# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(Shape_regularization_Tests)

cmake_minimum_required(VERSION 3.1...3.15)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)
  message(STATUS "Found CGAL")

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Find OSQP library and headers.
  find_package(OSQP REQUIRED)
  if(OSQP_FOUND)
    message(STATUS "Found OSQP")

    set(targets
    test_neighbor_query
    test_parallel_groups
    test_collinear_groups
    test_orthogonal_groups
    test_0_segments
    test_1_segment
    test_2_segments
    test_3_segments
    test_4_segments
    test_100_segments_angles
    test_100_segments_offsets
    test_cgal_solver
    test_plane_regularization
    test_closed_contour_1
    test_closed_contour_2
    test_open_contour
    test_equal_contours
    test_longest_direction
    test_multiple_directions_1
    test_multiple_directions_2
    test_user_defined_directions)

    foreach(target ${targets})
      create_single_source_cgal_program("${target}.cpp")
      if(TARGET ${target})
        target_link_libraries(${target} PUBLIC ${OSQP_LIBRARIES})
      endif()
    endforeach()

  else()
    message(STATUS "NOTICE: OSQP was not found. OSQP tests won't be available.")
  endif()
else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
