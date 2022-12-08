#================================================================
# cmake utilities to build python-bindings for kernel component
#================================================================

# wrap - Set extra dependencies
# set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
#   ${SICONOS_SWIG_SOURCE_DIR}/KernelTypes.i
#   KernelRegistration.i
#   ${SICONOS_SWIG_SOURCE_DIR}/solverOptions.i
#   ${SICONOS_SWIG_SOURCE_DIR}/handleException.i
#   ${SICONOS_SWIG_SOURCE_DIR}/start.i
#   ${SICONOS_SWIG_SOURCE_DIR}/ignored_functions.i
#   ${SICONOS_SWIG_SOURCE_DIR}/end.i
#   ${SICONOS_SWIG_SOURCE_DIR}/numpy.i
#   ${SICONOS_SWIG_SOURCE_DIR}/pyRegister.i
#   ${SICONOS_SWIG_SOURCE_DIR}/graph.i)

include(swig_python_tools)


## FP : when serialization is ON, we need includes from io in kernel ...
# This has to be reviewed.
if(WITH_SERIALIZATION)
  list(APPEND SWIG_KERNEL_INCLUDES ${CMAKE_SOURCE_DIR}/io/src/serialization)
  
  if(WITH_GENERATION)
    list(APPEND SWIG_KERNEL_INCLUDES  ${CMAKE_BINARY_DIR}/io/)
  else()
    list(APPEND SWIG_KERNEL_INCLUDES  ${CMAKE_SOURCE_DIR}/io/src/generation)
  endif()
  # for f2c.h in SiconosFull.h ...
  list(APPEND SWIG_KERNEL_INCLUDES  ${CMAKE_SOURCE_DIR}/externals/tools/)
endif()


add_swig_sub_module(
  FILE kernel.i
  INCLUDES ${SWIG_KERNEL_INCLUDES}
  DEPS kernel numerics)


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


