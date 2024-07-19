# -- OpenCASCADE --
# Must :
# - check if OpenCASCADE is available on the system
# Will:
# - add OpenCASCADE to the build process
# - set SICONOS_HAS_OpenCASCADE var (distributed in siconosConfig.cmake)

find_package(OpenCASCADE  REQUIRED)

if(OpenCASCADE_FOUND)
  message(STATUS "Found OpenCASCADE version ${OpenCASCADE_VERSION}")
  message("    OpenCASCADE libraries : ${OpenCASCADE_LIBRARIES}.")
  message("    OpenCASCADE headers path : ${OpenCASCADE_INCLUDE_DIR}")
endif()

if(NOT OpenCASCADE_VERSION VERSION_GREATER 7.6)
  message(FATAL_ERROR "Uncompatible opencascade version. Minimal required version is 7.6")
endif()

set(SICONOS_HAS_OpenCASCADE TRUE CACHE INTERNAL "True if OpenCASCADE API has been found and is activated.)")



if(OpenCASCADE_FOUND)

  if(OpenCASCADE_WITH_VTK) # Required for 7.6 but not for 7.7, 7.8  ... 
    find_package(VTK  REQUIRED  COMPONENTS CommonCore RenderingOpenGL2 RenderingFreeType)
  endif()

  if(NOT TARGET OpenCASCADE::OpenCASCADE)
    add_library(OpenCASCADE::OpenCASCADE IMPORTED INTERFACE)
    set_property(TARGET OpenCASCADE::OpenCASCADE PROPERTY INTERFACE_LINK_LIBRARIES ${OpenCASCADE_LIBRARIES})
    if(OpenCASCADE_INCLUDE_DIR)
      set_target_properties(OpenCASCADE::OpenCASCADE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${OpenCASCADE_INCLUDE_DIR}")
    endif()
  endif()


endif()


