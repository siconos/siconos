# -- OpenCASCADE --
# Must :
# - check if OpenCASCADE is available on the system
# Will:
# - add OpenCASCADE to the build process
# - set SICONOS_HAS_OpenCASCADE var (distributed in siconosConfig.cmake)

find_package(OpenCASCADE 7.4 REQUIRED)

if(OpenCASCADE_FOUND)
  message(STATUS "Found OpenCASCADE version ${OpenCASCADE_VERSION}")
  message("    OpenCASCADE libraries : ${OpenCASCADE_LIBRARIES}.")
  message("    OpenCASCADE headers path : ${OpenCASCADE_INCLUDE_DIR}")
endif()

set(SICONOS_HAS_OpenCASCADE TRUE CACHE INTERNAL "True if OpenCASCADE API has been found and is activated.)")

if(OpenCASCADE_FOUND)
  if(NOT TARGET OpenCASCADE::OpenCASCADE)
    add_library(OpenCASCADE::OpenCASCADE IMPORTED INTERFACE)
    set_property(TARGET OpenCASCADE::OpenCASCADE PROPERTY INTERFACE_LINK_LIBRARIES ${OpenCASCADE_LIBRARIES})
    if(OpenCASCADE_INCLUDE_DIR)
      set_target_properties(OpenCASCADE::OpenCASCADE PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${OpenCASCADE_INCLUDE_DIR}")
    endif()
  endif()
endif()
