# Find the Siconos Mechanics includes and libraries.
# The following variables are set if Mechanics is found.  If Mechanics is not
# found, SiconosMechanics_FOUND is set to false.
#  SiconosMechanics_FOUND        - True when the SiconosMechanics include directory is found.
#  SiconosMechanics_INCLUDE_DIRS - the path to where the Siconos Mechanics include files are.
#  SiconosMechanics_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosMechanics_LIBRARIES    - The libraries to link against Siconos Mechanics

# One may want to use a specific Mechanics Library by setting
# SiconosMechanics_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosMechanics)

IF(SiconosMechanics_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosMechanics_FOUND SiconosMechanics PATHS "${SiconosMechanics_LIBRARY_DIRECTORY}")
ELSE(SiconosMechanics_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosMechanics_FOUND SiconosMechanics)
ENDIF(SiconosMechanics_LIBRARY_DIRECTORY)

IF(SiconosMechanics_FOUND)
  GET_FILENAME_COMPONENT(SiconosMechanics_LIBRARY_DIRS ${SiconosMechanics_FOUND} PATH)
  SET(SiconosMechanics_LIBRARIES ${SiconosMechanics_FOUND})
  GET_FILENAME_COMPONENT(SiconosMechanics_LIBRARY_DIRS_DIR ${SiconosMechanics_LIBRARY_DIRS} PATH)
  SET(SiconosMechanics_INCLUDE_DIRS ${SiconosMechanics_LIBRARY_DIRS_DIR}/../include/Siconos/Mechanics)
ELSE(SiconosMechanics_FOUND)
  IF(SiconosMechanics_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Mechanics library not found. Please specify library location in SiconosMechanics_LIBRARY_DIRECTORY")
  ENDIF(SiconosMechanics_FIND_REQUIRED)
ENDIF(SiconosMechanics_FOUND)
