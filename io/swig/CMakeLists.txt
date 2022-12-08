#================================================================
# cmake utilities to build python-bindings for io component
#================================================================


if(HAVE_SICONOS_MECHANICS)

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

    # -- collect python files for the siconos/io module --
    collect_files(VAR pyfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR}/io EXTS py)
    if(NOT SICONOS_HAS_OCE)
      list(REMOVE_ITEM pyfiles io/occtools.py)
    endif()

    foreach(pyfile IN LISTS pyfiles)
      configure_file(${pyfile} ${SICONOS_SWIG_ROOT_DIR}/${pyfile} @ONLY)
    endforeach()

    # -- collect python files used to create scripts --

    collect_files(VAR scriptfiles DIRS ${CMAKE_CURRENT_SOURCE_DIR}/io/scripts EXTS py)

    foreach(pyfile IN LISTS scriptfiles)
      get_filename_component(short_name ${pyfile} NAME_WE)
      configure_file(${pyfile} ${SICONOS_SWIG_ROOT_DIR}/io/${short_name} @ONLY)
      install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/${short_name} DESTINATION bin)
    endforeach()

    # automatic manpages using help2man
    # See if help2man is available
    find_program(HELP2MAN help2man)
    if (HELP2MAN)
      # all other components must be built first so that tools don't fail to run
      set(_OTHER_COMPONENTS ${COMPONENTS})
      list(REMOVE_ITEM _OTHER_COMPONENTS ${COMPONENT})
      list(APPEND _OTHER_COMPONENTS  pyio_base)
      set(MANPAGES_FOR_TOOLS)
      include(tools4tests)

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
          DEPENDS ${SICONOS_SWIG_ROOT_DIR}/io/${TOOL} pyio_base 
          WORKING_DIRECTORY ${SICONOS_SWIG_ROOT_DIR}/..
          COMMENT "Generating man/${TOOL}.1 using help2man"
          )
        install(FILES ${CMAKE_BINARY_DIR}/man/${TOOL}.1 DESTINATION share/man/man1)
        list(APPEND MANPAGES_FOR_TOOLS ${CMAKE_BINARY_DIR}/man/${TOOL}.1)
      endmacro()

      
      if (HAVE_SICONOS_MECHANICS)
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
  if(INSTALL_PYTHON_SYMLINKS)
    message("Making symlinks for io Python modules")
    execute_process(COMMAND sh -c "mkdir -p '${SICONOS_SWIG_ROOT_DIR}/io'")
    configure_file(io/pprocess.py   ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess @ONLY)
  else()
    configure_file(io/pprocess.py   ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess @ONLY)
    install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/siconos_pprocess DESTINATION bin RENAME siconos_pprocess)
  endif()
endif()

include(swig_python_tools)


set(SWIG_IO_INCLUDES ${CMAKE_SOURCE_DIR}/kernel/swig)
set(SWIG_IO_DEPS numerics kernel io)
if(HAVE_SICONOS_MECHANICS)
  set(SWIG_IO_COMPILE_DEFINITIONS WITH_MECHANICS)
  list(APPEND SWIG_IO_INCLUDES  ${CMAKE_SOURCE_DIR}/io/src/mechanics)
  list(APPEND SWIG_IO_DEPS mechanics)
endif()



add_swig_sub_module(
  FILE io/io_base.i
  DEPS ${SWIG_IO_DEPS}
  INCLUDES ${SWIG_IO_INCLUDES}
  COMPILE_DEFINITIONS ${SWIG_IO_COMPILE_DEFINITIONS})

if(WITH_TESTING)
  if(NOT WITH_SERIALIZATION)
    list(APPEND python_excluded_tests tests/test_serialization.py)
  endif()
  build_python_tests(
    EXCLUDE "${python_excluded_tests}"
    )
endif()


