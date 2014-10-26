# Find the Siconos Control includes and libraries.
# The following variables are set if Control is found.  If Control is not
# found, SiconosControl_FOUND is set to false.
#  SiconosControl_FOUND        - True when the SiconosControl include directory is found.
#  SiconosControl_INCLUDE_DIRS - the path to where the Siconos Control include files are.
#  SiconosControl_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosControl_LIBRARIES    - The libraries to link against Siconos Control

# One may want to use a specific Control Library by setting
# SiconosControl_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosControl)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosControl_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for Control library in ${SiconosControl_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosControl_LIBRARY SiconosControl PATHS "${SiconosControl_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosControl_LIBRARY)
    MESSAGE(STATUS "Found : ${SiconosControl_LIBRARY}")
  ENDIF(SiconosControl_LIBRARY)
ELSE(SiconosControl_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosControl_LIBRARY SiconosControl ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
ENDIF(SiconosControl_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosControl
  REQUIRED_VARS SiconosControl_LIBRARY)

SET(SiconosControl_FOUND ${SICONOSCONTROL_FOUND})

IF(SiconosControl_LIBRARY)
  GET_FILENAME_COMPONENT(SiconosControl_LIBRARY_DIRS ${SiconosControl_LIBRARY} PATH)
  SET(SiconosControl_LIBRARIES ${SiconosControl_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosControl_LIBRARY_DIRS_DIR ${SiconosControl_LIBRARY_DIRS} PATH)

  FIND_PATH(SiconosControl_INCLUDE_DIRS ControlManager.hpp
    HINTS ${SiconosControl_LIBRARY_DIRS_DIR} ${SiconosControl_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES include/Siconos/Control)

  IF(NOT SiconosControl_INCLUDE_DIRS)
    IF(SiconosControl_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos Control headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
    ENDIF(SiconosControl_FIND_REQUIRED)
  ENDIF(NOT SiconosControl_INCLUDE_DIRS)

ELSE(SiconosControl_LIBRARY)
  IF(SiconosControl_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Control library not found. Please specify library location in SiconosControl_LIBRARY_DIRECTORY")
  ENDIF(SiconosControl_FIND_REQUIRED)
ENDIF(SiconosControl_LIBRARY)
