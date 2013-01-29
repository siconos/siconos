# Find the Siconos ACEF includes and libraries.
# The following variables are set if ACEF is found.  If ACEF is not
# found, SiconosACEF_FOUND is set to false.
#  SiconosACEF_FOUND        - True when the SiconosACEF include directory is found.
#  SiconosACEF_INCLUDE_DIRS - the path to where the Siconos ACEF include files are.
#  SiconosACEF_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosACEF_LIBRARIES    - The libraries to link against Siconos ACEF

# One may want to use a specific ACEF Library by setting
# SiconosACEF_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosACEF)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosACEF_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosACEF_LIBRARY SiconosACEF PATHS "${SiconosACEF_LIBRARY_DIRECTORY}")
ELSE(SiconosACEF_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosACEF_LIBRARY SiconosACEF)
ENDIF(SiconosACEF_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosACEF
  REQUIRED_VARS SiconosACEF_LIBRARY)

SET(SiconosACEF_FOUND ${SICONOSACEF_FOUND})

IF(SiconosACEF_LIBRARY)
  GET_FILENAME_COMPONENT(SiconosACEF_LIBRARY_DIRS ${SiconosACEF_LIBRARY} PATH)
  SET(SiconosACEF_LIBRARIES ${SiconosACEF_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosACEF_LIBRARY_DIRS_DIR ${SiconosACEF_LIBRARY_DIRS} PATH)
  SET(SiconosACEF_INCLUDE_DIRS ${SiconosACEF_LIBRARY_DIRS_DIR}/include/Siconos/ACEF)
ELSE(SiconosACEF_LIBRARY)
  IF(SiconosACEF_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos ACEF library not found. Please specify library location in SiconosACEF_LIBRARY_DIRECTORY")
  ENDIF(SiconosACEF_FIND_REQUIRED)
ENDIF(SiconosACEF_LIBRARY)
