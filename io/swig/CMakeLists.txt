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
  message(STATUS "input   io/SimpleGui.py ")
  message(STATUS "output ${SICONOS_SWIG_ROOT_DIR}/io/SimpleGui.py")


  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/io/swig/io/ io/img/*.*)
  message(STATUS "------------------------------ ${rendererFiles}")
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/io/swig/io/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      message(STATUS "          ")
      #message(STATUS "Configuring file ${rendererFile}")
      
      #message(STATUS "input   ${srcRendererPath}")
      #message(STATUS "output  ${SICONOS_SWIG_ROOT_DIR}/io/${rendererFile}")
      # configure_file(
      #   ${srcRendererPath}
      # 	${SICONOS_SWIG_ROOT_DIR}/io/${rendererFile}
      #   @ONLY)
      message(STATUS "Installing file ${rendererFile}")      
      message(STATUS "input   ${srcRendererPath}")
      message(STATUS "output  bin/img   ")    
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/img)
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)

  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/io/swig/io/ io/threeJS_libraries/*.*)
  message(STATUS "------------------------------ ${rendererFiles}")
  foreach(rendererFile ${rendererFiles})
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/io/swig/io/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      message(STATUS "          ")
      #message(STATUS "Configuring file ${rendererFile}")
      
      #message(STATUS "input   ${srcRendererPath}")
      #message(STATUS "output  ${SICONOS_SWIG_ROOT_DIR}/io/${rendererFile}")
      # configure_file(
      #   ${srcRendererPath}
      # 	${SICONOS_SWIG_ROOT_DIR}/io/${rendererFile}
      #   @ONLY)
      message(STATUS "Installing file ${rendererFile}")      
      message(STATUS "input   ${srcRendererPath}")
      message(STATUS "output  bin/threeJS_libraries")    
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/threeJS_libraries  )
    endif(NOT IS_DIRECTORY ${srcRendererPath})
  endforeach(rendererFile)  

  #configure_file(io/img ${SICONOS_SWIG_ROOT_DIR}/io/img @ONLY)
  #configure_file(io/threeJS_libraries ${SICONOS_SWIG_ROOT_DIR}/io/threeJS_libraries  @ONLY)

  install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/vview.py
    DESTINATION bin RENAME vview)
  install(PROGRAMS ${SICONOS_SWIG_ROOT_DIR}/io/pprocess.py
    DESTINATION bin RENAME pprocess)
  install(PROGRAMS io/renderer.py
    DESTINATION bin RENAME renderer)
endif()
include(swig_python_tools)
swig_module_setup(${COMPONENT}_PYTHON_MODULES)
build_python_tests()


