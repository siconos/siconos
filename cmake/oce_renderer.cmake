# json and python files only.
# This part is dedicated to the installation of oce renderer files (no build, only conf and installation).
if(HAVE_SICONOS_MECHANICS AND WITH_RENDERER)
 
  # This file is not installed as a program, even though it is one.
  # Remove the comment if you know how it should be installed --xhub
  #  configure_file(io/SimpleGui.py ${SICONOS_SWIG_ROOT_DIR}/io/SimpleGui.py @ONLY) 

  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals/renderer  ${CMAKE_SOURCE_DIR}/externals/renderer/img/*.*)
  foreach(rendererFile IN LISTS rendererFiles)
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/renderer/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/renderer/img)
    endif()
  endforeach()

  file(GLOB rendererFiles RELATIVE ${CMAKE_SOURCE_DIR}/externals/renderer/  ${CMAKE_SOURCE_DIR}/externals/renderer/threeJS_libraries/*.*)
  foreach(rendererFile IN LISTS rendererFiles)
    set(srcRendererPath ${CMAKE_SOURCE_DIR}/externals/renderer/${rendererFile})
    if(NOT IS_DIRECTORY ${srcRendererPath})     
      install(FILES  ${srcRendererPath} DESTINATION share/siconos/renderer/threeJS_libraries  )
    endif()
  endforeach()

  if(INSTALL_PYTHON_SYMLINKS)
    message("Setting up symlink install targets for externals Python executables")
    install(CODE "execute_process(COMMAND sh -c \"test -e '${CMAKE_INSTALL_PREFIX}/bin/siconos_renderer' || ln -vs '${CMAKE_CURRENT_SOURCE_DIR}/renderer/renderer.py' '${CMAKE_INSTALL_PREFIX}/bin/siconos_renderer'\")")
  else()
    install(PROGRAMS renderer/renderer.py
      DESTINATION bin RENAME siconos_renderer)
  endif()

endif()
