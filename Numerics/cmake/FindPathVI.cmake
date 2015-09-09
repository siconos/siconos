# Find the Path Ferris includes and libraries.
# The following variables are set if Ferris is found.  If Ferris is not
# found, PathVI_FOUND is set to false.
#  PathVI_FOUND        - True when the PathVI include directory is found.
#  PathVI_INCLUDE_DIRS - the path to where the Path Ferris include files are.
#  PathVI_LIBRARY_DIRS - The path to where the Path library files are.
#  PathVI_LIBRARIES    - The libraries to link against Path Ferris

# One may want to use a specific Ferris Library by setting
# PathVI_LIBRARY_DIRECTORY before FIND_PACKAGE(PathVI)
INCLUDE(FindPackageHandleStandardArgs)

IF(NOT PathVI_VERSION)
  SET(PathVI_VERSION "")
ENDIF(NOT PathVI_VERSION)

IF(NOT PathVI_LIBRARY_DIRECTORY AND GAMS_FOUND)
  SET(PathVI_LIBRARY_DIRECTORY ${GAMS_DIR})
ENDIF(NOT PathVI_LIBRARY_DIRECTORY AND GAMS_FOUND)

IF(PathVI_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathVI_LIBRARY pathvi${PathVI_VERSION} PATHS "${PathVI_LIBRARY_DIRECTORY}")
ELSE(PathVI_LIBRARY_DIRECTORY)
  FIND_LIBRARY(PathVI_LIBRARY pathvi${PathVI_VERSION})
ENDIF(PathVI_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(PathVI
  REQUIRED_VARS PathVI_LIBRARY)

SET(PathVI_FOUND ${PATHVI_FOUND})

IF(PathVI_LIBRARY)
  #  GET_FILENAME_COMPONENT(PathVI_LIBRARY_DIRS ${PathVI_LIBRARY} PATH)
  # GET_FILENAME_COMPONENT(PathVI_LIBRARIES ${PathVI_LIBRARY} NAME)
  #  GET_FILENAME_COMPONENT(PathVI_LIBRARY_DIRS_DIR ${PathVI_LIBRARY_DIRS} PATH)
  #  SET(PathVI_INCLUDE_DIRS ${PathVI_LIBRARY_DIRS_DIR}/include)
  SET(HAVE_PATHVI TRUE)
ELSE(PathVI_LIBRARY)
  IF(PathVI_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required path${PathVI_VERSION} library not found. Please specify library location in PathVI_LIBRARY_DIRECTORY")
  ENDIF(PathVI_FIND_REQUIRED)

ENDIF(PathVI_LIBRARY)
