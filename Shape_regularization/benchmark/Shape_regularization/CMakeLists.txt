# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

project(Shape_regularization_Benchmarks)

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
    benchmark_contours
    benchmark_qp_segments)

    foreach(target ${targets})
      create_single_source_cgal_program("${target}.cpp")
      if(TARGET ${target})
      target_link_libraries(${target} PUBLIC ${OSQP_LIBRARIES})
      endif()
    endforeach()

  else()
    message(STATUS "NOTICE: OSQP was not found. OSQP benchmarks won't be available.")
  endif()
else()
  message(WARNING
  "This program requires the CGAL library, and will not be compiled.")
endif()
