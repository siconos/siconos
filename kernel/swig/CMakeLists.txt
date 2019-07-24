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

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)

set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
if(NOT WITH_XML)
  list(APPEND ${COMPONENT}_python_excluded_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_xml.py)
endif()

if(NOT HAS_FORTRAN)
  list(APPEND ${COMPONENT}_python_excluded_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_matrix_exp.py)
endif()

build_python_tests()

configure_file(fromXml.py.in ${SICONOS_SWIG_ROOT_DIR}/fromXml.py @ONLY)
