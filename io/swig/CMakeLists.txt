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

  if(INSTALL_PYTHON_SYMLINKS)
    message("Making symlinks for io Python modules")
    execute_process(COMMAND sh -c "mkdir -p '${SICONOS_SWIG_ROOT_DIR}/io'")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/mechanics_io.py' '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py'")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/FrictionContactTrace.py' ${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py")

    message("Setting up symlink install targets for io Python executables")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_vview' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/vview.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_vview' \")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_vexport' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/vexport.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_vexport'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_pprocess' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/pprocess.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_pprocess'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_info' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/info.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_info'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_filter' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/filter.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_filter'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_run' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/run.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_run'\")")
  else()
    configure_file(io/mechanics_io.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_io.py @ONLY)
    configure_file(io/FrictionContactTrace.py ${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py @ONLY)
    install(PROGRAMS io/vview.py    DESTINATION bin RENAME siconos_vview)
    install(PROGRAMS io/vexport.py  DESTINATION bin RENAME siconos_vexport)
    install(PROGRAMS io/pprocess.py DESTINATION bin RENAME siconos_pprocess)
    install(PROGRAMS io/info.py     DESTINATION bin RENAME siconos_info)
    install(PROGRAMS io/filter.py   DESTINATION bin RENAME siconos_filter)
    install(PROGRAMS io/run.py      DESTINATION bin RENAME siconos_run)
  endif()

endif()
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

# right now the only tests we have require serialization
IF(WITH_SERIALIZATION)
  build_python_tests()
ENDIF(WITH_SERIALIZATION)


