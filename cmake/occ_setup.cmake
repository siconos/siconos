# -- OpenCASCADE --
# Must :
# - check if oce is available on the system
# Will:
# - add occ to the build process
# - set SICONOS_HAS_OpenCASCADE var (distributed in siconosConfig.cmake)

if(WITH_OpenCASCADE)
  # Set the list of toolkits modules which are requested.
  # Do we really need to set this explicitely??
  set(OpenCASCADE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")
    set(OpenCASCADE_TOOLKITS "TKernel"  "TKMath")
  # Explore system to find oce and required toolkits --> Looks for OpenCASCADEConfig.cmake and OpenCASCADEConfigVersion.cmake which
  # must have been installed with oce package.
  # It defines the following variables (according to OpenCASCADE doc)
  #  OpenCASCADE_INCLUDE_DIRS - include directory for OpenCASCADE
  #  OpenCASCADE_LIBRARIES    - all libraries to link against (warning, may be slower than just specify the used libs)
  #  OpenCASCADE_ALL_FOUND    - set to TRUE if all requested COMPONENTS are specified (see below), false otherwise
  #  OpenCASCADE_MISSING_TOOLKITS - when OpenCASCADE_ALL_FOUND is FALSE, contains a list of missing toolkits
  #  OpenCASCADE_ALL_BUILT_MODULES - the list of source directories compiled (mostly useful when running swig to generate wrappers)
  #if(OpenCASCADE_ALL_FOUND)
  #  message(STATUS "OpenCASCADE found.")
  #  message("    OpenCASCADE libraries : ${OpenCASCADE_LIBRARIES}.")
  #  message("    OpenCASCADE headers path : ${OpenCASCADE_INCLUDE_DIRS}.")
  #  message(STATUS "Found OpenCASCADE version ${OpenCASCADE_VERSION}")
  #else()
  #  message(FATAL_ERROR "OpenCASCADE detection failed due to missing toolkit(s): ${OpenCASCADE_MISSING_TOOLKITS}")
  #endif()
  
  #find_package(OpenCASCADE REQUIRED COMPONENTS ${OpenCASCADE_TOOLKITS})
  find_package(OpenCASCADE)
  
  

  message("############# OpenCASCADE found." ${OpenCASCADE_FOUND}, ${OpenCASCADE_INSTALL_PREFIX})
  message(STATUS "OpenCASCADE found.")
  message("    OpenCASCADE libraries : ${OpenCASCADE_LIBRARIES}.")
  message("    OpenCASCADE headers path : ${OpenCASCADE_INCLUDE_DIR}.")
  message(STATUS "Found OpenCASCADE version ${OpenCASCADE_VERSION}")
 
		 

  # if(OpenCASCADE_ALL_FOUND)
  #   message(STATUS "OpenCASCADE found.")
  #   message("    OpenCASCADE libraries : ${OpenCASCADE_LIBRARIES}.")
  #   message("    OpenCASCADE headers path : ${OpenCASCADE_INCLUDE_DIRS}.")
  #   message(STATUS "Found OpenCASCADE version ${OpenCASCADE_VERSION}")
  # else()
  #   message(FATAL_ERROR "OpenCASCADE detection failed due to missing toolkit(s): ${OpenCASCADE_MISSING_TOOLKITS}")
  # endif()


  set(SICONOS_HAS_OpenCASCADE TRUE CACHE INTERNAL "True if OpenCASCADE API has been found and is activated.)")

  # For versions of OpenCASCADE older than 0.18 AND on some specific systems (namely Debian),
  # some toolkits are detected (through OpenCASCADEConfig.cmake) but the
  # libraries are not available. We must not link them with the component.
  #if(OpenCASCADE_VERSION VERSION_LESS 0.18)
  set(UNNEEDED_OpenCASCADE_TOOLKITS "DRAWEXE" "TKDraw" "TKTopTest" "TKViewerTest" "TKXSDRAW" "TKDCAF" "TKXDEDRAW" "TKTObjDRAW" "TKQADraw")
  foreach(_T IN LISTS UNNEEDED_OpenCASCADE_TOOLKITS)
    list(REMOVE_ITEM OpenCASCADE_LIBRARIES  ${_T})
  endforeach()
  #endif()

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
endif()
