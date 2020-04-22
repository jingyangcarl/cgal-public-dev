# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(Shape_regularization_Examples)

cmake_minimum_required(VERSION 2.8.10)
set(CMAKE_CXX_STANDARD 11)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  # Use Eigen.
  find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)

  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})

    # Find OSQP library and headers.
    find_package(osqp REQUIRED)
    if(osqp_FOUND)
      message(STATUS "Found OSQP")

      set(targets 
      regularize_framework 
      regularize_15_segments 
      regularize_100_segments_angles
      regularize_100_segments_offsets
      regularize_simple
      regularize_real_data_2)
      
      foreach(target ${targets})
        create_single_source_cgal_program("${target}.cpp")
        if(TARGET ${target})
          target_link_libraries(${target} PUBLIC osqp::osqp)
        endif()
      endforeach()

    else()
      message(STATUS "NOTICE: OSQP was not found. OSQP examples won't be available.")
    endif()

    create_single_source_cgal_program(
      "regularize_planes.cpp")

  else()
    message(WARNING 
    "This program requires the Eigen3 library, and will not be compiled.")
  endif()

  create_single_source_cgal_program(
    "regularize_closed_contour.cpp")
  create_single_source_cgal_program(
    "regularize_open_contour.cpp")

else()
  message(WARNING 
  "This program requires the CGAL library, and will not be compiled.")
endif()
