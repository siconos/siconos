# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

include(tools4tests)

if(WITH_TESTING)
  add_custom_target(kernel-tests echo "Start kernel tests")

  # ---- Siconos Algebra tests ----
  begin_tests(src/algebra/test)

  new_test(
    NAME testSiconosAlgebra
    SOURCES  SiconosMatrixTest.cpp SiconosVectorTest.cpp BlockVectorTest.cpp ${SIMPLE_TEST_MAIN}
    DEPS "numerics;CPPUNIT::CPPUNIT;externals"
    )
  new_test(
    NAME testSiconosMemory
    SOURCES SiconosMemoryTest.cpp ${SIMPLE_TEST_MAIN}
    DEPS "numerics;CPPUNIT::CPPUNIT"
    )


  # ---- Siconos tools tests ----
  begin_tests(src/utils/test DEPS "CPPUNIT::CPPUNIT")
  new_test(SOURCES SiconosGraphTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES  SiconosPropertiesTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES  QuaternionTest.cpp ${SIMPLE_TEST_MAIN} DEPS numerics)

  # # ---- Modeling tools ---
  begin_tests(src/modeling/test DEPS "numerics;CPPUNIT::CPPUNIT")
  new_test(SOURCES FirstOrderNonLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearDSTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearTIRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderLinearRTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES FirstOrderType1RTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianScleronomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianRheonomousRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianCompliantRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianCompliantLinearTIRTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianSparseRelationsTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianSparseDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianSparseLinearTIDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES LagrangianLinearTIDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES NewtonEulerDSTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES NewtonEulerRelationsTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES RelationCompatibilityTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES NonSmoothDynamicalSystemTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES SiconosVisitorTest.cpp  ${SIMPLE_TEST_MAIN})
  new_test(SOURCES BoundaryConditionTest.cpp  ${SIMPLE_TEST_MAIN})

  # # ---- Simulation tools ---
  begin_tests(src/simulation/test DEPS "externals;numerics;CPPUNIT::CPPUNIT")
  new_test(SOURCES LsodarTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES EulerMoreauTest.cpp ${SIMPLE_TEST_MAIN})
  new_test(SOURCES OSNSPTest.cpp ${SIMPLE_TEST_MAIN})
  # Non-regression tests of the time-stepping schemes (refs: <test name>.ref in
  # src/simulation/test, regenerated with SICONOS_REGENERATE_TEST_REFS=1)
  new_test(SOURCES TimeSteppingSchemesTest.cpp ${SIMPLE_TEST_MAIN}
    COMPILE_DEFINITIONS "SICONOS_TEST_REF_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/src/simulation/test\"")

  # Optional validation of the time-stepping schemes against hots-python (Python port of
  # the Scilab code Hots, Acary CMAME 2013), see src/simulation/test/validate_with_hots.py.
  # hots-python is downloaded at configure time, unless HOTS_PYTHON_DIR points to a clone.
  option(WITH_HOTS_VALIDATION
    "Compare the time-stepping schemes with hots-python (downloaded at configure time)" OFF)
  if(WITH_HOTS_VALIDATION)
    set(HOTS_PYTHON_GIT_REPOSITORY
      "https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/hots-python.git"
      CACHE STRING "git repository of hots-python (an ssh URL works too)")
    set(HOTS_PYTHON_GIT_TAG "main" CACHE STRING "branch, tag or commit of hots-python")
    set(HOTS_PYTHON_DIR "" CACHE PATH
      "existing hots-python source directory (no download when set)")
    if(NOT HOTS_PYTHON_DIR)
      include(FetchContent)
      message(STATUS "hots-python: ${HOTS_PYTHON_GIT_REPOSITORY} (${HOTS_PYTHON_GIT_TAG})")
      FetchContent_Declare(hots_python
        GIT_REPOSITORY ${HOTS_PYTHON_GIT_REPOSITORY}
        GIT_TAG ${HOTS_PYTHON_GIT_TAG}
        GIT_SHALLOW TRUE
        UPDATE_DISCONNECTED TRUE # no git fetch at each configure
        SOURCE_SUBDIR no_cmake_project) # python only: nothing to add_subdirectory
      FetchContent_MakeAvailable(hots_python)
      set(_hots_dir ${hots_python_SOURCE_DIR})
    else()
      set(_hots_dir ${HOTS_PYTHON_DIR})
    endif()
    if(NOT EXISTS ${_hots_dir}/hots/__init__.py)
      message(FATAL_ERROR "WITH_HOTS_VALIDATION: no hots package in ${_hots_dir}")
    endif()
    if(NOT Python_EXECUTABLE)
      find_package(Python COMPONENTS Interpreter REQUIRED)
    endif()
    message(STATUS "hots-python validation with ${_hots_dir}")

    # 1. run the scheme tests and dump the full trajectories (<test>.dat) in the test
    #    directory (where the .ref files are copied) ...
    set(_hots_work ${CMAKE_CURRENT_BINARY_DIR}/src/simulation/test)
    add_test(NAME TimeSteppingSchemesTest_dump
      COMMAND $<TARGET_FILE:TimeSteppingSchemesTest>
      WORKING_DIRECTORY ${_hots_work})
    set_tests_properties(TimeSteppingSchemesTest_dump PROPERTIES
      ENVIRONMENT "SICONOS_TEST_DUMP_TRAJECTORIES=1"
      FIXTURES_SETUP hots_trajectories)
    # 2. ... and integrate the same systems with hots-python
    add_test(NAME TimeSteppingSchemes_vs_hots
      COMMAND ${Python_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation/test/validate_with_hots.py
        --hots ${_hots_dir} --dir ${_hots_work}
      WORKING_DIRECTORY ${_hots_work})
    set_tests_properties(TimeSteppingSchemes_vs_hots PROPERTIES
      FIXTURES_REQUIRED hots_trajectories)
  endif()

  if(LPSOLVE_FOUND)
    # HAS_EXTREME_POINT_ALGO must be added to compile def if lpsolve has been found.
    list(APPEND compile_defs HAS_EXTREME_POINT_ALGO)
  endif()

  new_test(SOURCES testAVI.cpp ${SIMPLE_TEST_MAIN} DEPS LAPACK::LAPACK COMPILE_DEFINITIONS ${compile_defs})
  if(HAS_FORTRAN)
    new_test(SOURCES ZOHTest.cpp ${SIMPLE_TEST_MAIN} DEPS LAPACK::LAPACK)
  endif()

 endif()
