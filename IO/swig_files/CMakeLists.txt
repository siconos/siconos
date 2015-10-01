#================================================================
# cmake utilities to build python-bindings for io component
#================================================================

# python_bindings - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "./${COMPONENT}")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)

set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()

if(HAVE_SICONOS_MECHANICS)
  configure_file(mechanics_io.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py @ONLY)
  configure_file(vview.py ${SICONOS_SWIG_ROOT_DIR}/io/vview.py @ONLY)

  install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/vview.py
    DESTINATION bin RENAME vview)
  
endif()

