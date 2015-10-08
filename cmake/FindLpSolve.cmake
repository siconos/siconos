# -----------------------------------------
# Search for lpsolve header and library
# --> lp_lib.h
# --> liblpsolve
#
# User hint : LpSolve_DIR --> try cmake -DLpSolve_DIR=path-to-lpsolve-install
#
# header is searched for in 'standard' path + optional suffix lpsolve, LpSolve_DIR/include
# LpSolve_DIR/lpsolve
#
# Please note that we are using lpsolve 5.5
# -----------------------------------------
include(LibFindMacros)
#
# Header
# first, search in user-defined LpSolve_DIR
if(LpSolve_DIR)
  find_path(LpSolve_INCLUDE_DIR
    NAMES lp_lib.h
    PATHS ${LpSolve_DIR}
    PATH_SUFFIXES lpsolve include
    NO_DEFAULT_PATH)
endif()
# else, standard path
find_path(
  LpSolve_INCLUDE_DIR
  NAMES lp_lib.h
  PATH_SUFFIXES lpsolve include
  )

IF (NOT LpSolve_INCLUDE_DIR)
  MESSAGE(STATUS "Cannot find LPSOLVE headers")
ELSE()
  if(LpSolve_DIR)
    find_library(LpSolve_LIBRARY lpsolve55
      PATHS ${Lp_Solve_DIR}
      PATH_SUFFIXES lib lib/${CMAKE_LIBRARY_ARCHITECTURE} 
      )
  endif()
  # debian nonsense: see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=503314
  if(EXISTS "/etc/debian_version")
    FIND_LIBRARY(LpSolve_LIBRARY lpsolve55 PATH_SUFFIXES "lp_solve")
    IF(LpSolve_LIBRARY AND "${LpSolve_LIBRARY}" MATCHES ".a$" AND NOT I_WANT_STATIC_LPSOLVE)
      SET(LpSolve_LIBRARY FALSE)
    ENDIF()
  ELSE()
    FIND_LIBRARY(LpSolve_LIBRARY
      lpsolve55
      PATHS ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
      )
  ENDIF()
  
  IF (LpSolve_LIBRARY)
    GET_FILENAME_COMPONENT(LpSolve_LIBRARY_DIRS ${LpSolve_LIBRARY} PATH)
    GET_FILENAME_COMPONENT(LpSolve_LIBRARIES ${LpSolve_LIBRARY} REALPATH)
  #ELSE(LpSolve_LIBRARY)
  #  MESSAGE("Cannot find LPSOLVE libraries!")
  endif()
  IF (NOT LpSolve_FIND_QUIETLY)
    MESSAGE(STATUS "Found LpSolve: ${LpSolve_LIBRARY}")
    MESSAGE(STATUS "LpSolve_LIBRARIES: ${LpSolve_LIBRARIES}")
  ENDIF (NOT LpSolve_FIND_QUIETLY)
  
ENDIF (NOT LpSolve_INCLUDE_DIR)

set(LpSolve_PROCESS_LIBS LpSolve_LIBRARIES)
set(LpSolve_PROCESS_INCLUDES LpSolve_INCLUDE_DIR)
libfind_process(LpSolve)
