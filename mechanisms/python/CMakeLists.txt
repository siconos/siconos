#==================================================================
# cmake utilities to build python-bindings for mechanisms component
#==================================================================

## PB11 wrap for mechanisms NOT YET IMPLEMENTED - TO BE DONE

# -- install the script --
set(EXAMPLE_PATH ./) # @EXAMPLE_PATH@ is required in siconos_mechanisms script
configure_file(siconos_mechanisms.py.in ${CMAKE_BINARY_DIR}/scripts/siconos_mechanisms.py)
configure_file(mbtbDefaultOptions.py ${CMAKE_BINARY_DIR}/python/siconos/mechanisms/mbtbDefaultOptions.py)
install(PROGRAMS ${CMAKE_BINARY_DIR}/scripts/siconos_mechanisms.py DESTINATION bin RENAME siconos_mechanisms)


# --- Tests ---
if(WITH_TESTING)
  # -- prepare siconos_mechanisms tests --
  set(EXAMPLE_PATH ${CMAKE_BINARY_DIR}/python/tests/mechanisms) # @EXAMPLE_PATH@ is required in siconos_mechanisms script
  # Each file is copied into siconos/tests.
  configure_file(siconos_mechanisms.py.in  ${EXAMPLE_PATH}/siconos_mechanisms.py)

  # data (CAD and body) taken from slider crank example in mechanisms component
  set(path_to_test_sources ${CMAKE_SOURCE_DIR}/mechanisms/python/tests/slider_crank)
  configure_file(${path_to_test_sources}/bodydef.py ${EXAMPLE_PATH}/bodydef.py COPYONLY)
  configure_file(${path_to_test_sources}/mbtbLocalOptions.py ${EXAMPLE_PATH} COPYONLY)
  
  build_python_tests(
    # DEPS Boost::boost # for quaternions
    )
endif()
