#================================================================
# cmake utilities to build python-bindings for io component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "io/io_base")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)

set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")
set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/kernel/swig)

if(HAVE_SICONOS_MECHANICS)
  configure_file(io/mechanics_io.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py @ONLY)
  configure_file(io/vview.py ${SICONOS_SWIG_ROOT_DIR}/io/vview.py @ONLY)
  configure_file(io/pprocess.py ${SICONOS_SWIG_ROOT_DIR}/io/pprocess.py @ONLY)
  configure_file(io/SimpleGui.py ${SICONOS_SWIG_ROOT_DIR}/io/SimpleGui.py @ONLY)

  install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/vview.py
    DESTINATION bin RENAME vview)
  install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/pprocess.py
    DESTINATION bin RENAME pprocess)
endif()
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()


