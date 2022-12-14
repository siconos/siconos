#================================================================
# cmake utilities to build python-bindings for control component
#================================================================

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES
  "control/controller.i"
  "control/observer.i"
  "control/sensor.i"
  "control/simulation.i"
)

include(swig_python_tools)
foreach(tgt IN LISTS ${COMPONENT}_PYTHON_MODULES)

  add_swig_sub_module(
    FILE ${tgt}
    DEPS control
    INCLUDES ${CMAKE_SOURCE_DIR}/kernel/swig
    )
endforeach()

if(WITH_TESTING)
  if(NOT HAS_FORTRAN)
    list(APPEND python_excluded_tests tests/test_smc.py)
  endif()
  build_python_tests(
    # DEPS ${COMPONENT} # some plugins in tests need to be linked with kernel
    EXCLUDE ${python_excluded_tests})
endif()
