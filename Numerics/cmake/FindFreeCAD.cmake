# Find FreeCAD includes and libraries.
# The following variables are set if FreeCAD is found.  If FreeCAD is not
# found, FreeCAD_FOUND is set to false.
#  FreeCAD_FOUND        - True when the FreeCAD include directory is found.
#  FreeCAD_LIBRARY_DIRS - The path to where the Siconos library files are.
#  FreeCAD_LIBRARIES    - The libraries to link against Siconos FreeCAD

# One may want to use a specific FreeCAD Library by setting
# FreeCAD_LIBRARY_DIRECTORY before FIND_PACKAGE(FreeCAD)
INCLUDE(FindPackageHandleStandardArgs)

SET(_NEEDED_LIBRARIES FreeCAD ${FreeCAD_FIND_COMPONENTS})

FOREACH(_COMPONENT ${_NEEDED_LIBRARIES})
  IF(FreeCAD_LIBRARY_DIRECTORY)
    FIND_LIBRARY(FreeCAD_${_COMPONENT}_LIBRARY ${_COMPONENT}${CMAKE_SHARED_LIBRARY_SUFFIX} PATHS "${FreeCAD_LIBRARY_DIRECTORY}")
  ELSE(FreeCAD_LIBRARY_DIRECTORY)
    FIND_LIBRARY(FreeCAD_${_COMPONENT}_LIBRARY ${_COMPONENT}${CMAKE_SHARED_LIBRARY_SUFFIX} PATHS /usr/lib/freecad/lib /usr/local/lib/freecad/lib /opt/lib/freecad/lib)
  ENDIF()
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(FreeCAD
    REQUIRED_VARS FreeCAD_${_COMPONENT}_LIBRARY)
ENDFOREACH()



FOREACH(_COMPONENT ${_NEEDED_LIBRARIES})
  IF(NOT FreeCAD_${_COMPONENT}_LIBRARY)
    IF(FreeCAD_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required FreeCAD ${_COMPONENT}${CMAKE_SHARED_LIBRARY_SUFFIX} library not found. Please specify library location in FreeCAD_LIBRARY_DIRECTORY")
    ENDIF()
  ENDIF()
  LIST(APPEND FreeCAD_LIBRARIES ${FreeCAD_${_COMPONENT}_LIBRARY})
ENDFOREACH()

MESSAGE(STATUS "FreeCAD libraries : ${FreeCAD_LIBRARIES}")
SET(FreeCAD_FOUND TRUE)
