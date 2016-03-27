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
  configure_file(io/vview.py ${CMAKE_BINARY_DIR}/io/vview.py @ONLY)
  configure_file(io/pprocess.py ${CMAKE_BINARY_DIR}/io/pprocess.py @ONLY)
  # This files is not installed as a program, even though it is one.
  # Remove the comment if you know how it should be installed --xhub
  #  configure_file(io/SimpleGui.py ${SICONOS_SWIG_ROOT_DIR}/io/SimpleGui.py @ONLY) 


  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals  ${CMAKE_SOURCE_DIR}/externals/img/*.*)
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/img)
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)

  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals  ${CMAKE_SOURCE_DIR}/externals/threeJS_libraries/*.*)
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/threeJS_libraries  )
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)  

  install(PROGRAMS ${CMAKE_BINARY_DIR}/io/vview.py
    DESTINATION bin RENAME siconos_vview)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/io/pprocess.py
    DESTINATION bin RENAME siconos_pprocess)
  install(PROGRAMS io/renderer.py
    DESTINATION bin RENAME siconos_renderer)
endif()
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)

# right now the only tests we have require serialization
IF(WITH_SERIALIZATION)
  build_python_tests()
ENDIF(WITH_SERIALIZATION)


