#================================================================
# cmake utilities to build python-bindings for mechanics component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)

# DEFS
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "mechanics/bodies"
  "mechanics/joints"
  "mechanics/contact_detection/base")

if(WITH_OCC)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "mechanics/occ")
endif()

if(WITH_BULLET)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "mechanics/contact_detection/bullet")
  list(APPEND ${COMPONENT}_SWIG_DEFS "-DWITH_BULLET")
  if(BULLET_USE_DOUBLE_PRECISION)
    list(APPEND ${COMPONENT}_SWIG_DEFS "-DBT_USE_DOUBLE_PRECISION")
  endif()
endif()

if(WITH_MECHANISMS)
  list(APPEND ${COMPONENT}_PYTHON_MODULES
    "mechanics/mechanisms/cadmbtb"
    "mechanics/mechanisms/mbtb")
endif()

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/kernel/swig)

configure_file(mechanics/contact_detection/tools.py ${SICONOS_SWIG_ROOT_DIR}/mechanics/contact_detection/tools.py)
configure_file(mechanics/contact_detection/convexhull.py  ${SICONOS_SWIG_ROOT_DIR}/mechanics/contact_detection/convexhull.py )
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()


# siconos mechanisms script
if(WITH_MECHANISMS)
  # we need siconos-mechanism for build and install tree with different 'PYTHONPATH' for each of them.
  set(CURRENT_PYTHON_PATH ${CMAKE_BINARY_DIR}/wrap)
  configure_file(mechanics/mechanisms/siconos-mechanisms.py.in ${CMAKE_BINARY_DIR}/scripts_for_build/siconos_mechanisms.py)
  set(CURRENT_PYTHON_PATH ${SICONOS_PYTHON_INSTALL_DIR})
  configure_file(mechanics/mechanisms/siconos-mechanisms.py.in ${CMAKE_BINARY_DIR}/scripts/siconos_mechanisms.py)
  configure_file(mechanics/mechanisms/mbtbDefaultOptions.py ${SICONOS_SWIG_ROOT_DIR}/mechanics/mechanisms/mbtbDefaultOptions.py)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/scripts/siconos_mechanisms.py DESTINATION bin RENAME siconos_mechanisms)
endif()
