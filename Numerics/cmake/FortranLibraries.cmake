# cache result to avoid re-running
#
IF(NOT HAVE_FORTRAN_LIBRARIES)
# search for libraries required to link with fortran
# placed in the variable
#  FORTRAN_LIBRARIES
# the location is in
#  FORTRAN_COMPILER_LIB_DIRECTORY

# Need Fortran support for this
ENABLE_LANGUAGE(Fortran)

# these two can be used by the user to guide the process

GET_FILENAME_COMPONENT(fortran_comp_bin_dir ${CMAKE_Fortran_COMPILER} PATH)
GET_FILENAME_COMPONENT(fortran_comp_dir ${fortran_comp_bin_dir}       PATH)

SET(FORTRAN_COMPILER_DIRECTORY ${fortran_comp_dir} CACHE PATH "location of fortran compiler top directory")
SET(FORTRAN_COMPILER_LIB_DIRECTORY ${fortran_comp_dir}/lib CACHE PATH "location of fortran compiler lib directory")

SET(KNOWN_FORTRAN_LIBRARIES
   ""
   "gfortran"
   "fio\;f90math\;f77math"
   "afio\;af90math\;af77math\;amisc"
   "g2c"
)

# Define the list "options" of all possible schemes that we want to consider
# Get its length and initialize the counter "iopt" to zero
LIST(LENGTH KNOWN_FORTRAN_LIBRARIES imax)
SET(iopt 0)

# Try to link against fortran libraries
WHILE(${iopt} LESS ${imax})

  # Get the current list entry (current scheme)
  LIST(GET KNOWN_FORTRAN_LIBRARIES ${iopt} fc_libraries_t)
  
  STRING(REGEX REPLACE ";" "\\\;" fc_libraries "${fc_libraries_t}")
  
  # Create a simple Fortran library source
  FILE(WRITE ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB/ftest.f
    "        FUNCTION fortransub()\n"
    "        print *, 'hello world'\n"
    "        print *, ' value = ', atan(-1.0)\n"
    "        RETURN\n"
    "        END\n"
    )

  # a call to the Fortran library from C
  FILE(WRITE ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB/testF77libs.c
    "        int main(int argc, char ** argv) { fortransub_();}\n"
    )
  
  # let's find it in some potential libraries directories
  MESSAGE(STATUS "Seeking fortran library directory")
  
  FILE(GLOB _LIBDIRS_MAYBE 
    /opt /opt/* /opt/lib/*
    /usr/local/* /usr/local/lib/*
    /usr/* /usr/lib/*)
  
  FOREACH(_F ${_LIBDIRS_MAYBE})
    IF(IS_DIRECTORY ${_F})
      LIST(APPEND _LIBDIRS ${_F})
    ENDIF(IS_DIRECTORY ${_F})
  ENDFOREACH(_F ${_LIBDIRS_MAYBE})

  LIST(APPEND _LIBDIRS "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp")
  
  #string -> list
  SET(_libs ${fc_libraries})

  FIND_LIBRARY(_LFIND NAMES ${_libs} PATHS ${_LIBDIRS} PATH_SUFFIXES lib)
      
  IF(_LFIND)

    MESSAGE(STATUS "Found ${_LFIND}, trying to link mixed C/fortran stuff...")
    GET_FILENAME_COMPONENT(_LFINDDIR ${_LFIND} PATH)
    
    # Create a CMakeLists.txt file which will generate the "flib" library
    FILE(WRITE ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB/CMakeLists.txt
      "PROJECT(FortranTest C Fortran)\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n"
      "LINK_DIRECTORIES(${_LFINDDIR})\n"
      "ADD_LIBRARY(flib ftest.f)\n"
      "ADD_EXECUTABLE(testF77libs testF77libs.c)\n"
      "TARGET_LINK_LIBRARIES(testF77libs flib ${fc_libraries})\n")
    

    # Use TRY_COMPILE to make the target "flib"
    TRY_COMPILE(FORT_LIBS_WORK
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB
      flib
      OUTPUT_VARIABLE MY_OUTPUT)
    
    IF(FORT_LIBS_WORK)
      MESSAGE(STATUS "Ok, mixed C/Fortran can be linked with ${fc_libraries}")
      SET(FORTRAN_COMPILER_LIB_DIRECTORY ${_LFINDDIR})
    ENDIF(FORT_LIBS_WORK)
  ENDIF(_LFIND)
    
  IF(FORT_LIBS_WORK)
    # the goal is to set this to something useful
    SET(FORTRAN_LIBRARIES ${fc_libraries} CACHE STRING "fortran libraries required to link with f90 generated code" FORCE)
    MESSAGE(STATUS "Using FORTRAN_LIBRARIES ${FORTRAN_LIBRARIES}")
    MESSAGE(STATUS "Using FORTRAN_COMPILER_LIB_DIRECTORY ${FORTRAN_COMPILER_LIB_DIRECTORY}")
    
    SET(HAVE_FORTRAN_LIBRARIES TRUE CACHE INTERNAL "successfully found fortran libraries")
    SET(iopt ${imax})
    LINK_DIRECTORIES(${FORTRAN_COMPILER_LIB_DIRECTORY})
  ELSE(FORT_LIBS_WORK)
    MATH(EXPR iopt ${iopt}+1)
  ENDIF(FORT_LIBS_WORK)
ENDWHILE(${iopt} LESS ${imax})
  
IF(NOT HAVE_FORTRAN_LIBRARIES)
  MESSAGE(FATAL_ERROR "Could not find valid fortran link libraries. Try setting FORTRAN_COMPILER_DIRECTORY or FORTRAN_LIBRARIES")
ENDIF(NOT HAVE_FORTRAN_LIBRARIES)

ENDIF(NOT HAVE_FORTRAN_LIBRARIES)
