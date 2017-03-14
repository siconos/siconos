# -----------------------------------------
# Search for pathvi library
# --> libpathvi
#
# User hint : PathVI_DIR --> try cmake -DPathVI_DIR=path-to-pathvi-install
#
# header is searched for in 'standard' path + optional suffix , PathVI_DIR/include
# PathVI_DIR
#
# -----------------------------------------
include(LibFindMacros)

# ----------------------------------------
# First, search in (optionnaly) user-defined PathVI_DIR
# Try : 
# PathVI_DIR/
# PathVI_DIR/include
# ----------------------------------------
IF(NOT PathVI_VERSION)
  SET(PathVI_VERSION "")
ENDIF(NOT PathVI_VERSION)

IF(NOT PathVI_DIR AND GAMSCAPI_FOUND)
  SET(PathVI_DIR ${GAMS_DIR})
ENDIF(NOT PathVI_DIR AND GAMSCAPI_FOUND)

# ----------------------------------------
# Library == libpath${PathVI_VERSION}
# First, search in (optionnaly) user-defined PathVI_DIR
# Try:
# - PathVI_DIR
# - PathVI_DIR/lib
# - PathVI_DIR/lib/CMAKE_LIBRARY_ARCHITECTURE (for example on Debian like system : lib/x86_64-linux-gnu
if(PathVI_DIR)
find_library(PathVI_LIBRARY pathvi${PathVI_VERSION}
  PATHS ${PathVI_DIR}
  PATH_SUFFIXES lib lib/${CMAKE_LIBRARY_ARCHITECTURE} 
  )
else(PathVI_DIR)
find_library(PathVI_LIBRARY pathvi${PathVI_VERSION}
  )
endif(PathVI_DIR)
# If not found, try standard path.

# Set PathVI_LIBRARY_DIRS and PathVI_LIBRARIES for libfindprocess
# see Using LibFindMacros : https://cmake.org/Wiki/CMake:How_To_Find_Libraries
IF (PathVI_LIBRARY)
  GET_FILENAME_COMPONENT(PathVI_LIBRARY_DIRS ${PathVI_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(PathVI_LIBRARIES ${PathVI_LIBRARY} REALPATH)
  SET(HAVE_PATHVI TRUE)
endif()
IF (NOT PathVI_FIND_QUIETLY)
  MESSAGE(STATUS "Found PathVI: ${PathVI_LIBRARY}")
  MESSAGE(STATUS "PathVI_LIBRARIES: ${PathVI_LIBRARIES}")
ENDIF (NOT PathVI_FIND_QUIETLY)

# Final check :
set(PathVI_PROCESS_LIBS PathVI_LIBRARIES)
libfind_process(PathVI)
