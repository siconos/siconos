#================================================================
# cmake utilities to build python-bindings for kernel component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${SICONOS_SWIG_SRC_DIRS}/KernelTypes.i
  KernelRegistration.i
  ${SICONOS_SWIG_SRC_DIRS}/solverOptions.i
  ${SICONOS_SWIG_SRC_DIRS}/handleException.i
  ${SICONOS_SWIG_SRC_DIRS}/start.i
  ${SICONOS_SWIG_SRC_DIRS}/ignored_functions.i
  ${SICONOS_SWIG_SRC_DIRS}/end.i
  ${SICONOS_SWIG_SRC_DIRS}/numpy.i
  ${SICONOS_SWIG_SRC_DIRS}/pyRegister.i
  ${SICONOS_SWIG_SRC_DIRS}/graph.i)

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

target_link_libraries(_kernel numerics)

# --- Tests ---
if(WITH_${COMPONENT}_TESTING)
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
