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
if(NOT WITH_CXX)
  list(REMOVE_ITEM ${COMPONENT}_SWIG_DEFS "-dirvtable")
endif()

if(WITH_FCLIB)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "./fclib")
  # if(NOT FCLib_FCLIB_HEADER_ONLY)
  #   list(APPEND CMAKE_SWIG_FLAGS "-DFCLIB_NOT_HEADER_ONLY")
  # endif()
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
    set_property(SOURCE fclib.i PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)
  else() #USE_TARGET_INCLUDE property is not available in old cmake versions.
    find_package(FCLIB 3.0.0 CONFIG REQUIRED)
    list(APPEND CMAKE_SWIG_FLAGS "-I${FCLIB_INCLUDE_DIR}")
  endif()
else()
  list(APPEND python_excluded_tests tests/test_fclib.py)
endif()

if(NOT NO_RUNTIME_BUILD_DEP)
  # -- Scipy --
  # It is required for sparse matrices in numerics
  find_python_module(scipy  VERSION 0.15) # for sparse

  # -- mpi4py -- 
  if(WITH_MPI)
    find_python_module(mpi4py REQUIRED INCLUDES)
    set(WITH_MPI4PY TRUE)
    set(numerics_SWIG_DEFS "-I${mpi4py_INCLUDE_DIR};-DWITH_MPI4PY")

  endif()
endif()
include(swig_python_tools)

swig_module_setup(${COMPONENT}_PYTHON_MODULES)

# -- mpi4py -- 
if(WITH_MPI)
  target_include_directories(_numerics PRIVATE ${MPI_C_INCLUDE_DIRS})
  target_link_libraries(_numerics ${MPI_C_LIBRARIES})
  find_python_module(mpi4py REQUIRED INCLUDES)
  target_include_directories(_numerics PRIVATE ${mpi4py_INCLUDE_DIR})
  set(WITH_MPI4PY TRUE)
  target_compile_definitions(_numerics PRIVATE WITH_MPI4PY)
  #  "-I${mpi4py_INCLUDE_DIR};-DWITH_MPI4PY")
endif()


if(WITH_${COMPONENT}_TESTING)
  if(WITH_MUMPS AND WITH_MPI4PY)
    message(STATUS "add test python mumps + mpi4py")
  else()
    list(APPEND python_excluded_tests tests/test_NM_MUMPS.py)
  endif()

  include(tools4tests)
  build_python_tests(
    DEPS ${COMPONENT} # some plugins in tests need to be linked with numerics
    EXCLUDE ${python_excluded_tests})
endif()
