#================================================================
# cmake utilities to build python-bindings for numerics component
#================================================================

if(NOT NO_RUNTIME_BUILD_DEP)
  # -- Scipy --
  # It is required for sparse matrices in numerics
  find_python_module(scipy  REQUIRED)# for sparse

endif()

include(swig_python_tools)


list(APPEND deps $<BUILD_INTERFACE:externals> numerics)
if(WITH_SYSTEM_SUITESPARSE)
  list(APPEND SWIG_NUMERICS_DEPS $<BUILD_INTERFACE:SuiteSparse::CXSparse>)
endif()


# -- mpi4py -- 
if(WITH_MPI)
  list(APPEND SWIG_NUMERICS_INCLUDES ${MPI_C_INCLUDE_DIRS})
  list(APPEND SWIG_NUMERICS_DEPS ${MPI_C_LIBRARIES})
  
  find_python_module(mpi4py REQUIRED INCLUDES)
  set(WITH_MPI4PY TRUE)
  list(APPEND SWIG_NUMERICS_INCLUDES ${mpi4py_INCLUDE_DIR})
  list(APPEND SWIG_NUMERICS_COMPILE_DEFINITIONS WITH_MPI4PY)
endif()

if(WITH_SERIALIZATION)
  list(APPEND SWIG_NUMERICS_INCLUDES ${CMAKE_SOURCE_DIR}/io/src/serialization)
endif()


add_swig_sub_module(
  FILE numerics.i
  DEPS ${SWIG_NUMERICS_DEPS} numerics externals
  INCLUDES ${SWIG_NUMERICS_INCLUDES}
  COMPILE_DEFINITIONS ${SWIG_NUMERICS_COMPILE_DEFINITIONS})


if(WITH_FCLIB)
  add_swig_sub_module(
    FILE fclib.i
    DEPS FCLIB::fclib numerics # fclib swig depends on numerics because of start.i. This is quite illogical ... 
    )
else()
  list(APPEND python_excluded_tests tests/test_fclib.py)
endif()


if(WITH_TESTING)
  if(WITH_MUMPS AND WITH_MPI4PY)
    message(STATUS "add test python mumps + mpi4py")
  else()
    list(APPEND python_excluded_tests tests/test_NM_MUMPS.py)
  endif()
  
  include(tools4tests)
  if(WITH_SYSTEM_SUITESPARSE)
    set(suitesparse SuiteSparse::CXSparse)
  else()
    set(suitesparse)
  endif()

  build_python_tests(
    DEPS "${COMPONENT};${suitesparse}" # some plugins in tests need to be linked with numerics
    EXCLUDE ${python_excluded_tests}
    DEPS externals)
endif()
