# Find the Siconos Kernel includes and libraries.
# The following variables are set if Kernel is found.  If Kernel is not
# found, SiconosKernel_FOUND is set to false.
#  SiconosKernel_FOUND        - True when the SiconosKernel include directory is found.
#  SiconosKernel_INCLUDE_DIRS - the path to where the Siconos Kernel include files are.
#  SiconosKernel_LIBRARY_DIRS - The path to where the Siconos library files are.
#  SiconosKernel_LIBRARIES    - The libraries to link against Siconos Kernel

# One may want to use a specific Kernel Library by setting
# SiconosKernel_LIBRARY_DIRECTORY before FIND_PACKAGE(SiconosKernel)

IF(SiconosKernel_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosKernel_FOUND SiconosKernel PATHS "${SiconosKernel_LIBRARY_DIRECTORY}")
ELSE(SiconosKernel_LIBRARY_DIRECTORY)
  FIND_LIBRARY(SiconosKernel_FOUND SiconosKernel)
ENDIF(SiconosKernel_LIBRARY_DIRECTORY)

IF(SiconosKernel_FOUND)
  GET_FILENAME_COMPONENT(SiconosKernel_LIBRARY_DIRS ${SiconosKernel_FOUND} PATH)
  SET(SiconosKernel_LIBRARIES ${SiconosKernel_FOUND})
  GET_FILENAME_COMPONENT(SiconosKernel_LIBRARY_DIRS_DIR ${SiconosKernel_LIBRARY_DIRS} PATH)
  SET(SiconosKernel_INCLUDE_DIRS ${SiconosKernel_LIBRARY_DIRS_DIR}/include/Siconos/Kernel)
ELSE(SiconosKernel_FOUND)
  IF(SiconosKernel_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required Siconos Kernel library not found. Please specify library location in SiconosKernel_LIBRARY_DIRECTORY")
  ENDIF(SiconosKernel_FIND_REQUIRED)
ENDIF(SiconosKernel_FOUND)
