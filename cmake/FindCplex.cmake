# -----------------------------------------
# Search for cplex library
# --> libcplex
#
# User hint : Cplex_DIR --> try cmake -DCplex_DIR=path-to-cplex-install
#
# header is searched for in 'standard' path + optional suffix , Cplex_DIR/include
# Cplex_DIR
#
# -----------------------------------------
include(LibFindMacros)

# ----------------------------------------
# First, search in (optionnaly) user-defined Cplex_DIR
# Try : 
# Cplex_DIR/
# Cplex_DIR/include
# ----------------------------------------

# XXX we do not look for include now ... --xhub

IF(NOT Cplex_DIR AND GAMS_FOUND)
  SET(Cplex_DIR ${GAMS_DIR})
ENDIF(NOT Cplex_DIR AND GAMS_FOUND)

# ----------------------------------------
# Library == libpath${Cplex_VERSION}
# First, search in (optionnaly) user-defined Cplex_DIR
# Try:
# - Cplex_DIR
# - Cplex_DIR/lib
# - Cplex_DIR/lib/CMAKE_LIBRARY_ARCHITECTURE (for example on Debian like system : lib/x86_64-linux-gnu
if(Cplex_DIR)
find_library(Cplex_LIBRARY cplex
  PATHS ${Cplex_DIR}
  PATH_SUFFIXES lib lib/${CMAKE_LIBRARY_ARCHITECTURE} 
  )
else(Cplex_DIR)
find_library(Cplex_LIBRARY cplex
  )
endif(Cplex_DIR)
# If not found, try standard path.

# Set Cplex_LIBRARY_DIRS and Cplex_LIBRARIES for libfindprocess
# see Using LibFindMacros : https://cmake.org/Wiki/CMake:How_To_Find_Libraries
IF (Cplex_LIBRARY)
  GET_FILENAME_COMPONENT(Cplex_LIBRARY_DIRS ${Cplex_LIBRARY} PATH)
  GET_FILENAME_COMPONENT(Cplex_LIBRARIES ${Cplex_LIBRARY} REALPATH)
endif()
IF (NOT Cplex_FIND_QUIETLY)
  MESSAGE(STATUS "Found Cplex: ${Cplex_LIBRARY}")
  MESSAGE(STATUS "Cplex_LIBRARIES: ${Cplex_LIBRARIES}")
ENDIF (NOT Cplex_FIND_QUIETLY)

# Final check :
set(Cplex_PROCESS_LIBS Cplex_LIBRARIES)
libfind_process(Cplex)
