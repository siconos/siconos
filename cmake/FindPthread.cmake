# Find external lpthread includes and libraries.
# The following variables are set if liblpthread is found.  If liblpthread is not
# found, Pthread_FOUND is set to false.
#  Pthread_FOUND        - True when the Pthread include directory is found.
#  Pthread_INCLUDE_DIRS - the path to where the Pthread include files are.
#  Pthread_LIBRARY_DIRS - The path to where the Pthread library files are.
#  Pthread_LIBRARIES    - The libraries to link against Pthread

# One may want to use a specific Ferris Library by setting
# Pthread_LIBRARY_DIRECTORY before FIND_PACKAGE(Pthread)
INCLUDE(FindPackageHandleStandardArgs)

IF(Pthread_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Pthread_LIBRARY pthread PATHS "${Pthread_LIBRARY_DIRECTORY}")
ELSE(Pthread_LIBRARY_DIRECTORY)
  FIND_LIBRARY(Pthread_LIBRARY pthread)
ENDIF(Pthread_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pthread
  REQUIRED_VARS Pthread_LIBRARY)

SET(Pthread_FOUND ${PTHREAD_FOUND})

IF(Pthread_LIBRARY)
  GET_FILENAME_COMPONENT(Pthread_LIBRARY_DIRS ${Pthread_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(Pthread_LIBRARIES ${Pthread_LIBRARY} NAME)
  GET_FILENAME_COMPONENT(Pthread_LIBRARY_DIRS_DIR ${Pthread_LIBRARY_DIRS} PATH)
  # SET(Pthread_INCLUDE_DIRS ${Pthread_LIBRARY_DIRS_DIR}/include)  Note FP: failed on ubuntu because of suffix after lib ... and useless anyway.
ELSE(Pthread_LIBRARY)
  IF(Pthread_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Pthread library not found. Please specify library location in Pthread_LIBRARY_DIRECTORY")
  ENDIF(Pthread_FIND_REQUIRED)

ENDIF(Pthread_LIBRARY)
