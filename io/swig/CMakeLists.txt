#================================================================
# cmake utilities to build python-bindings for io component
#================================================================

# wrap - Set extra dependencies
set(SWIG_MODULE_${COMPONENT}_EXTRA_DEPS)

# List of modules to be built in the python bindings
set(${COMPONENT}_PYTHON_MODULES "io/io_base")

set(${COMPONENT}_SWIG_INCLUDE_DIRECTORIES)

set(${COMPONENT}_SWIG_DEFS "${SWIG_DEFS}")

list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/kernel/swig)
if(HAVE_SICONOS_MECHANICS)
  list(APPEND ${COMPONENT}_SWIG_DEFS "-DWITH_MECHANICS")

  if(INSTALL_PYTHON_SYMLINKS)
    message("Making symlinks for io Python modules")
    execute_process(COMMAND sh -c "mkdir -p '${SICONOS_SWIG_ROOT_DIR}/io'")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_hdf5.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/mechanics_hdf5.py' '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_hdf5.py'")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_run.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/mechanics_run.py' '${SICONOS_SWIG_ROOT_DIR}/io/mechanics_run.py'")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/FrictionContactTrace.py' ${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py")
    execute_process(COMMAND sh -c "test -e '${SICONOS_SWIG_ROOT_DIR}/io/vview.py' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/vview.py' ${SICONOS_SWIG_ROOT_DIR}/io/vview.py")

    message("Setting up symlink install targets for io Python executables")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_vview' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/siconos_vview.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_vview' \")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_vexport' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/siconos_vexport.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_vexport'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_export_raw_data' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/siconos_export_export_raw_data.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_export_raw_data'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_info' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/info.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_info'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_filter' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/filter.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_filter'\")")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_run' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/io/run.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_run'\")")
  else()

    configure_file(io/mechanics_hdf5.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_hdf5.py @ONLY)
    configure_file(io/mechanics_run.py ${SICONOS_SWIG_ROOT_DIR}/io/mechanics_run.py @ONLY)
    configure_file(io/FrictionContactTrace.py ${SICONOS_SWIG_ROOT_DIR}/io/FrictionContactTrace.py @ONLY)
    configure_file(io/vview.py      ${SICONOS_SWIG_ROOT_DIR}/io/vview.py @ONLY)

    configure_file(io/siconos_vview.py ${SICONOS_SWIG_ROOT_DIR}/io/siconos_vview @ONLY)
    configure_file(io/siconos_vexport.py ${SICONOS_SWIG_ROOT_DIR}/io/siconos_vexport @ONLY)
    configure_file(io/siconos_export_raw_data.py ${SICONOS_SWIG_ROOT_DIR}/io/siconos_export_raw_data @ONLY)
    configure_file(io/info.py       ${SICONOS_SWIG_ROOT_DIR}/io/siconos_info @ONLY)
    configure_file(io/filter.py     ${SICONOS_SWIG_ROOT_DIR}/io/siconos_filter @ONLY)
    configure_file(io/run.py        ${SICONOS_SWIG_ROOT_DIR}/io/siconos_run @ONLY)
    configure_file(io/compare.py    ${SICONOS_SWIG_ROOT_DIR}/io/siconos_compare @ONLY)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_vview    DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_vexport  DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_export_raw_data  DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_info     DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_filter   DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_run      DESTINATION bin)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_compare  DESTINATION bin)
    # automatic manpages using help2man

    # See if help2man is available
    find_program(HELP2MAN help2man)
    if (HELP2MAN)
      # all other components must be built first so that tools don't fail to run
      set(_OTHER_COMPONENTS ${COMPONENTS})
      list(REMOVE_ITEM _OTHER_COMPONENTS ${COMPONENT})

      # swig modules must be built before vview and vexport can generate manpages
      foreach(_SWIG_MODULE ${${COMPONENT}_PYTHON_MODULES})
        get_filename_component(_SWIG_MODULE_NAME ${_SWIG_MODULE} NAME)
        list(APPEND _OTHER_COMPONENTS ${_SWIG_MODULE_NAME}_swig_compilation _${_SWIG_MODULE_NAME})
      endforeach()

      set(MANPAGES_FOR_TOOLS)
      include(tools4tests)
      set_ldlibpath()
      macro (GEN_MANPAGE_FROM_HELP TOOL)
        # set PYTHONPATH to swig dir, which contains "siconos", so
        # that "import siconos.*" works
        add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/man/${TOOL}.1
          COMMAND env PYTHONPATH=${SICONOS_SWIG_ROOT_DIR}/.. ${LDLIBPATH}
            ${HELP2MAN} -s1 -N -o  ${CMAKE_BINARY_DIR}/man/${TOOL}.1
            --opt-include=${CMAKE_SOURCE_DIR}/man/${TOOL}.help2man
            --include=${CMAKE_SOURCE_DIR}/man/footer.help2man
            --no-discard-stderr
            -n \"`sed '/Description: /!d\;s///' ${SICONOS_SWIG_ROOT_DIR}/io/${TOOL}`\"
            ${SICONOS_SWIG_ROOT_DIR}/io/${TOOL}
          DEPENDS ${SICONOS_SWIG_ROOT_DIR}/io/${TOOL} ${_OTHER_COMPONENTS}
          WORKING_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/..
          COMMENT "Generating man/${TOOL}.1 using help2man"
          )
        install(FILES ${CMAKE_BINARY_DIR}/man/${TOOL}.1 DESTINATION share/man/man1)
        list(APPEND MANPAGES_FOR_TOOLS ${CMAKE_BINARY_DIR}/man/${TOOL}.1)
      endmacro()
      list(FIND COMPONENTS mechanics _HAVE_MECHANICS)
      if (NOT ${_HAVE_MECHANICS} EQUAL -1)
        GEN_MANPAGE_FROM_HELP(siconos_vview)
        GEN_MANPAGE_FROM_HELP(siconos_vexport)
        GEN_MANPAGE_FROM_HELP(siconos_export_raw_data)
        GEN_MANPAGE_FROM_HELP(siconos_info)
        GEN_MANPAGE_FROM_HELP(siconos_run)
        GEN_MANPAGE_FROM_HELP(siconos_filter)
        GEN_MANPAGE_FROM_HELP(siconos_compare)
      endif()

      add_custom_target(manpages ALL DEPENDS ${MANPAGES_FOR_TOOLS})
    endif()
  endif()

endif()

if(HAVE_SICONOS_MECHANISMS)
  list(APPEND ${COMPONENT}_SWIG_DEFS "-DWITH_MECHANISMS")
  if(${CMAKE_VERSION} VERSION_LESS "3.13")
    list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/mechanisms/src/")
  endif()
  if(INSTALL_PYTHON_SYMLINKS)
    message("Making symlinks for io Python modules")
    execute_process(COMMAND sh -c "mkdir -p '${SICONOS_SWIG_ROOT_DIR}/io'")
    configure_file(io/pprocess.py   ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess @ONLY)
  else()
    configure_file(io/pprocess.py   ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess @ONLY)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess DESTINATION bin RENAME siconos_pprocess)
  endif()
endif()

# TMP fix 
if(${CMAKE_VERSION} VERSION_LESS "3.13")
  # Notice that with last cmake versions, it's possible to propagate includes to swig in a better way.
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/numerics/src/tools")
  list(APPEND ${COMPONENT}_SWIG_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIR}")
endif()




include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

# right now the only tests we have require serialization
if(WITH_SERIALIZATION)
  if(WITH_${COMPONENT}_TESTING)
    build_python_tests(
      # DEPS ${COMPONENT} # some plugins in tests need to be linked with kernel
      )
  endif()
endif()


