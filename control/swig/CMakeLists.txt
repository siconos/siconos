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

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/kernel/swig)
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")
# TMP fix 
if(${CMAKE_VERSION} VERSION_LESS "3.13")
  # Notice that with last cmake versions, it's possible to propagate includes to swig in a better way.
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/numerics/src/tools/")
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}")
endif()

include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

if(WITH_TESTING)
  if(NOT HAS_FORTRAN)
    list(APPEND python_excluded_tests tests/test_smc.py)
  endif()
  build_python_tests(
    # DEPS ${COMPONENT} # some plugins in tests need to be linked with kernel
    EXCLUDE ${python_excluded_tests})
endif()
