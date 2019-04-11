#================================================================
# cmake utilities to build python-bindings for numerics component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${SICONOS_SWIG_SRC_DIRS}/solverOptions.i
  ${SICONOS_SWIG_SRC_DIRS}/start.i
  ${SICONOS_SWIG_SRC_DIRS}/end.i
  ${SICONOS_SWIG_SRC_DIRS}/numpy.i)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "./${COMPONENT}")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)
# numerics is not written in C++
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")
IF(NOT WITH_CXX)
  list(REMOVE_ITEM ${COMPONENT}_SWIG_DEFS "-dirvtable")
ENDIF(NOT WITH_CXX)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
#list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
#if("${isSystemDir}" STREQUAL "-1")
#  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#endif("${isSystemDir}" STREQUAL "-1")

if(WITH_FCLIB)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "./fclib")
  IF(NOT FCLib_FCLIB_HEADER_ONLY)
  list(APPEND CMAKE_SWIG_FLAGS "-DFCLIB_NOT_HEADER_ONLY")
  ENDIF()
else()
  list(APPEND ${COMPONENT}_python_excluded_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_fclib.py)
endif()

if(WITH_${COMPONENT}_TESTING)
  if(WITH_MUMPS AND WITH_MPI4PY)
    message(STATUS "add test python mumps + mpi4py")
  else()
    list(APPEND ${COMPONENT}_python_excluded_tests
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_NM_MUMPS.py)
  endif()
endif()

if(NOT NO_RUNTIME_BUILD_DEP)
  # -- Scipy --
  # It is required for sparse matrices in numerics
  find_python_module(scipy  VERSION 0.15) # for sparse
endif()
include(swig_python_tools)

swig_module_setup(${COMPONENT}_PYTHON_MODULES)

build_python_tests()
