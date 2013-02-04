# Find the Siconos Mechanics includes and libraries.
# The following variables are set if Mechanics is found.  If Mechanics is not
# found, SiconosMechanics_FOUND is set to false.
#  SiconosMechanics_FOUND        - True when the SiconosMechanics include directory is found.
#  SiconosMechanics_INCLUDE_DIRS - the path to where the Siconos Mechanics include files are.
#  SiconosMechanics_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosMechanics_LIBRARIES    - The libraries to link against Siconos Mechanics

# One may want to use a specific Mechanics Library by setting
# SiconosMechanics_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosMechanics)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosMechanics_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for Mechanics library in ${SiconosMechanics_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosMechanics_LIBRARY SiconosMechanics PATHS "${SiconosMechanics_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosMechanics_LIBRARY)
    MESSAGE(STATUS "Found : ${SiconosMechanics_LIBRARY}")
  ENDIF(SiconosMechanics_LIBRARY)
ELSE(SiconosMechanics_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosMechanics_LIBRARY SiconosMechanics)
ENDIF(SiconosMechanics_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosMechanics
  REQUIRED_VARS SiconosMechanics_LIBRARY)

SET(SiconosMechanics_FOUND ${SICONOSMECHANICS_FOUND})

IF(SiconosMechanics_LIBRARY)
  GET_FILENAME_COMPONENT(SiconosMechanics_LIBRARY_DIRS ${SiconosMechanics_LIBRARY} PATH)
  SET(SiconosMechanics_LIBRARIES ${SiconosMechanics_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosMechanics_LIBRARY_DIRS_DIR ${SiconosMechanics_LIBRARY_DIRS} PATH)

   FIND_PATH(SiconosMechanics_INCLUDE_DIRS SpaceFilter.hpp
    HINTS ${SiconosMechanics_LIBRARY_DIRS_DIR} ${SiconosMechanics_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES Siconos/Mechanics)

  IF(NOT SiconosMechanics_INCLUDE_DIRS)
    IF(SiconosMechanics_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos Mechanics headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
    ENDIF(SiconosMechanics_FIND_REQUIRED)
  ENDIF(NOT SiconosMechanics_INCLUDE_DIRS)

ELSE(SiconosMechanics_LIBRARY)
  IF(SiconosMechanics_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Mechanics library not found. Please specify library location in SiconosMechanics_LIBRARY_DIRECTORY")
  ENDIF(SiconosMechanics_FIND_REQUIRED)
ENDIF(SiconosMechanics_LIBRARY)
