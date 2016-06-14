#================================================================
# cmake utilities to build python-bindings for mechanics component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)

# DEFS
set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES
  "mechanics/joints"
  "mechanics/collision/base"
  "mechanics/collision/native"
  "mechanics/collision/bodies")

if(WITH_OCC)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "mechanics/occ")
endif()

if(WITH_BULLET)
  list(APPEND ${COMPONENT}_PYTHON_MODULES "mechanics/collision/bullet")
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

configure_file(mechanics/collision/tools.py ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/tools.py)
configure_file(mechanics/collision/convexhull.py  ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/convexhull.py )
configure_file(mechanics/collision/__init__.py  ${SICONOS_SWIG_ROOT_DIR}/mechanics/collision/__init__.py )
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()


# siconos mechanisms script
if(WITH_MECHANISMS)
  configure_file(mechanics/mechanisms/siconos-mechanisms.py.in ${CMAKE_BINARY_DIR}/scripts/siconos-mechanisms.py)
  configure_file(mechanics/mechanisms/mbtbDefaultOptions.py ${SICONOS_SWIG_ROOT_DIR}/mechanics/mechanisms/mbtbDefaultOptions.py)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/scripts/siconos-mechanisms.py DESTINATION bin RENAME siconos_mechanisms)
endif()
