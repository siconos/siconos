# -- OCE --
# Must :
# - check if oce is available on the system
# Will:
# - add occ to the build process
# - set SICONOS_HAS_OCE var (distributed in siconosConfig.cmake)

if(SICONOS_HAS_OCE)
  # Job has already been done ... return
  return()
endif()

if(WITH_OCE)
  # Set the list of toolkits modules which are requested.
  # Do we really need to set this explicitely??
  set(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")
  
  # Explore system to find oce and required toolkits --> Looks for OCEConfig.cmake and OCEConfigVersion.cmake which
  # must have been installed with oce package.
  # It defines the following variables (according to OCE doc)
  #  OCE_INCLUDE_DIRS - include directory for OCE
  #  OCE_LIBRARIES    - all libraries to link against (warning, may be slower than just specify the used libs)
  #  OCE_ALL_FOUND    - set to TRUE if all requested COMPONENTS are specified (see below), false otherwise
  #  OCE_MISSING_TOOLKITS - when OCE_ALL_FOUND is FALSE, contains a list of missing toolkits
  #  OCE_ALL_BUILT_MODULES - the list of source directories compiled (mostly useful when running swig to generate wrappers)
  find_package(OCE 0.16 REQUIRED COMPONENTS ${OCE_TOOLKITS})

  if(OCE_ALL_FOUND)
    message(STATUS "OCE found.")
    message("    OCE libraries : ${OCE_LIBRARIES}.")
    message("    OCE headers path : ${OCE_INCLUDE_DIRS}.")
    message(STATUS "Found OCE version ${OCE_VERSION}")
  else()
    message(FATAL_ERROR "OCE detection failed due to missing toolkit(s): ${OCE_MISSING_TOOLKITS}")
  endif()
  
  set(SICONOS_HAS_OCE TRUE)

  # For versions of OCE older than 0.18 AND on some specific systems (namely Debian),
  # some toolkits are detected (through OCEConfig.cmake) but the
  # libraries are not available. We must not link them with the component.
  #if(OCE_VERSION VERSION_LESS 0.18)
    set(UNNEEDED_OCE_TOOLKITS "DRAWEXE" "TKDraw" "TKTopTest" "TKViewerTest" "TKXSDRAW" "TKDCAF" "TKXDEDRAW" "TKTObjDRAW" "TKQADraw"   )
    foreach(_T ${UNNEEDED_OCE_TOOLKITS})
      list(REMOVE_ITEM OCE_LIBRARIES  ${_T})
    endforeach()
  #endif()
endif()
