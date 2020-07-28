# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Solver_interface_Examples)

find_package(CGAL QUIET)
if(CGAL_FOUND)

  # Use Eigen.
  find_package(Eigen3 3.1.0) # (requires 3.1.0 or greater)
  if(EIGEN3_FOUND)
    include(${EIGEN3_USE_FILE})
    create_single_source_cgal_program("singular_value_decomposition.cpp")
    create_single_source_cgal_program("sparse_solvers.cpp")
  endif()

  create_single_source_cgal_program("diagonalize_matrix.cpp")
  create_single_source_cgal_program("mixed_integer_program.cpp")

  # Find OSQP library and headers.
  find_package(OSQP REQUIRED)
  if(OSQP_FOUND)
    message(STATUS "Found OSQP")

    set(targets
    osqp_quadratic_program)

    foreach(target ${targets})
      create_single_source_cgal_program("${target}.cpp")
      if(TARGET ${target})
        target_link_libraries(${target} PUBLIC ${OSQP_LIBRARIES})
      endif()
    endforeach()
  else()
    message(STATUS "NOTICE: OSQP was not found. OSQP examples won't be available.")
  endif()

  find_package(SCIP QUIET)
  if(SCIP_FOUND)
    include_directories(BEFORE ${SCIP_INCLUDE_DIRS})
    target_link_libraries(mixed_integer_program PRIVATE ${SCIP_LIBRARIES})
    target_compile_definitions(mixed_integer_program PRIVATE -DCGAL_USE_SCIP)
    message("SCIP found and used.")
  else()

    find_package(GLPK QUIET)
    if(GLPK_FOUND)
      include_directories(BEFORE ${GLPK_INCLUDE_DIR})
      target_link_libraries(mixed_integer_program PRIVATE ${GLPK_LIBRARIES})
      target_compile_definitions(mixed_integer_program PRIVATE -DCGAL_USE_GLPK)
      message("GLPK found and used.")
    else()
      message(STATUS
        "NOTICE : This project requires either SCIP or GLPK, and will not be compiled. "
        "Please provide either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")
    endif()
  endif()
else()
  message(STATUS
    "NOTICE: This program requires the CGAL library, and will not be compiled.")
endif()
