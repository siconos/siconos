#================================================================
# cmake utilities to build python-bindings for kernel component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${SICONOS_SWIG_SOURCE_DIR}/KernelTypes.i
  KernelRegistration.i
  ${SICONOS_SWIG_SOURCE_DIR}/solverOptions.i
  ${SICONOS_SWIG_SOURCE_DIR}/handleException.i
  ${SICONOS_SWIG_SOURCE_DIR}/start.i
  ${SICONOS_SWIG_SOURCE_DIR}/ignored_functions.i
  ${SICONOS_SWIG_SOURCE_DIR}/end.i
  ${SICONOS_SWIG_SOURCE_DIR}/numpy.i
  ${SICONOS_SWIG_SOURCE_DIR}/pyRegister.i
  ${SICONOS_SWIG_SOURCE_DIR}/graph.i)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "./${COMPONENT}")

set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")
# TMP fix 
if(${CMAKE_VERSION} VERSION_LESS "3.13")
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/numerics/src/tools")
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}")

  # Kernel dirs are required as include for swig in next components.
  foreach(_dir IN LISTS ${COMPONENT}_DIRS)
    list(APPEND CMAKE_SWIG_FLAGS "-I${CMAKE_SOURCE_DIR}/${COMPONENT}/${_dir}")
  endforeach()
  list(REMOVE_DUPLICATES CMAKE_SWIG_FLAGS)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} CACHE INTERNAL "")
endif()


include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

target_link_libraries(_kernel PUBLIC numerics)
target_include_directories(_kernel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}) # addons.hpp ...

## FP : when serialization is ON, we need includes from io in kernel ...
# This has to be reviewed.
if(WITH_SERIALIZATION)
  target_include_directories(_kernel PRIVATE ${CMAKE_SOURCE_DIR}/io/src/serialization)
  if(WITH_GENERATION)
    target_include_directories(_kernel PRIVATE ${CMAKE_BINARY_DIR}/io/)
  else()
    target_include_directories(_kernel PRIVATE ${CMAKE_SOURCE_DIR}/io/src/generation)
  endif()
  # for f2c.h in SiconosFull.h ...
  target_include_directories(_kernel PRIVATE ${CMAKE_SOURCE_DIR}/externals/tools/)
endif()

# --- Tests ---
if(WITH_TESTING)
  if(NOT WITH_XML)
    list(APPEND python_excluded_tests tests/test_xml.py)
  endif()
  
  if(NOT HAS_FORTRAN)
    list(APPEND python_excluded_tests tests/test_matrix_exp.py)
  endif()
  
  find_python_module(lxml REQUIRED)
  build_python_tests(
    EXCLUDE ${python_excluded_tests})
endif()

configure_file(fromXml.py.in ${SICONOS_SWIG_ROOT_DIR}/fromXml.py @ONLY)
