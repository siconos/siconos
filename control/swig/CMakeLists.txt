#================================================================
# cmake utilities to build python-bindings for control component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS
  ${CMAKE_SOURCE_DIR}/kernel/swig/kernel.i
  ControlBase.i
  ControlTypemaps.i
  )

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES
  "${COMPONENT}/controller"
  "${COMPONENT}/observer"
  "${COMPONENT}/sensor"
  "${COMPONENT}/simulation"
)

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/kernel/swig)
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
if(NOT WITH_FORTRAN)
  list(APPEND ${COMPONENT}_python_excluded_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_smc.py)
endif()
build_python_tests()
