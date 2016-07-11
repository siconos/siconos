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
  list(APPEND ${COMPONENT}_SWIG_DEFS "-DWITH_MECHANICS")
  configure_file(io/mechanics_io.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py @ONLY)
  configure_file(io/FrictionContactTrace.py ${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py @ONLY)
  
  configure_file(io/vview.py ${CMAKE_BINARY_DIR}/io/vview.py @ONLY)
  configure_file(io/vexport.py ${CMAKE_BINARY_DIR}/io/vexport.py @ONLY)
  configure_file(io/pprocess.py ${CMAKE_BINARY_DIR}/io/pprocess.py @ONLY)


  install(PROGRAMS ${CMAKE_BINARY_DIR}/io/vview.py
    DESTINATION bin RENAME siconos_vview)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/io/vexport.py
    DESTINATION bin RENAME siconos_vexport)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/io/pprocess.py
    DESTINATION bin RENAME siconos_pprocess)

  
endif()
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

# right now the only tests we have require serialization
IF(WITH_SERIALIZATION)
  build_python_tests()
ENDIF(WITH_SERIALIZATION)


