# Find the Siconos Kernel includes and libraries.
# The following variables are set if Kernel is found.  If Kernel is not
# found, SiconosKernel_FOUND is set to false.
#  SiconosKernel_FOUND        - True when the SiconosKernel include directory is found.
#  SiconosKernel_INCLUDE_DIRS - the path to where the Siconos Kernel include files are.
#  SiconosKernel_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosKernel_LIBRARIES    - The libraries to link against Siconos Kernel

# One may want to use a specific Kernel Library by setting
# SiconosKernel_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosKernel)
INCLUDE(FindPackageHandleStandardArgs)

IF(SiconosKernel_LIBRARY_DIRECTORY)
  MESSAGE(STATUS "Looking for Kernel library in ${SiconosKernel_LIBRARY_DIRECTORY}")
  FIND_LIBRARY(SiconosKernel_LIBRARY SiconosKernel PATHS "${SiconosKernel_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
  IF(SiconosKernel_LIBRARY)
    MESSAGE(STATUS "Found : ${SiconosKernel_LIBRARY}")
  ENDIF(SiconosKernel_LIBRARY)
ELSE(SiconosKernel_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosKernel_LIBRARY SiconosKernel  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
ENDIF(SiconosKernel_LIBRARY_DIRECTORY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(SiconosKernel
  REQUIRED_VARS SiconosKernel_LIBRARY)

SET(SiconosKernel_FOUND ${SICONOSKERNEL_FOUND})

IF(SiconosKernel_LIBRARY)
  SET(SiconosKernel_LIBRARIES ${SiconosKernel_LIBRARY})
  GET_FILENAME_COMPONENT(SiconosKernel_LIBRARY_DIRS ${SiconosKernel_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(SiconosKernel_LIBRARY_DIRS_DIR ${SiconosKernel_LIBRARY_DIRS} PATH)
  GET_FILENAME_COMPONENT(SiconosKernel_LIBRARY_DIRS_DIR_DIR ${SiconosKernel_LIBRARY_DIRS_DIR} PATH)

  FIND_PATH(SiconosKernel_INCLUDE_DIRS SiconosKernel.hpp
    HINTS ${SiconosKernel_LIBRARY_DIRS_DIR} ${SiconosKernel_LIBRARY_DIRS_DIR_DIR} 
    ENV PATH
    PATH_SUFFIXES Siconos/Kernel)
  
  IF(NOT SiconosKernel_INCLUDE_DIRS)
    IF(SiconosKernel_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos Kernel headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
    ENDIF(SiconosKernel_FIND_REQUIRED)
  ENDIF(NOT SiconosKernel_INCLUDE_DIRS)

ELSE(SiconosKernel_LIBRARY)
  IF(SiconosKernel_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Kernel library not found. Please specify library location in SiconosKernel_LIBRARY_DIRECTORY")
  ENDIF(SiconosKernel_FIND_REQUIRED)
ENDIF(SiconosKernel_LIBRARY)
