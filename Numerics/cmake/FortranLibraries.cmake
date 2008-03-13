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

SET(KNOWN_FORTRAN_LIBRARIES "gfortran\;gfortranbegin"
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
    /usr/local/*/lib/gcc/*/*
    /usr/lib/gcc/*/*
    /usr/local/lib/gcc/*/*
    /usr/* /usr/lib/*)
  
  FOREACH(_F ${_LIBDIRS_MAYBE})
    IF(IS_DIRECTORY ${_F})
      LIST(APPEND _LIBDIRS ${_F})
    ENDIF(IS_DIRECTORY ${_F})
  ENDFOREACH(_F ${_LIBDIRS_MAYBE})

  LIST(APPEND _LIBDIRS "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp")
  
  #string -> list
  SET(_libs ${fc_libraries})
  SET(_LFIND)
  FOREACH(_l ${_libs})
    MESSAGE(STATUS "Searching for the library ${_l}")
    FIND_LIBRARY(_LFIND_${_l} NAMES ${_l} PATHS ${_LIBDIRS} PATH_SUFFIXES lib)
    IF(_LFIND_${_l})
      MESSAGE(STATUS "Found ${_LFIND_${_l}}")
      SET(_LFIND TRUE)
    ENDIF(_LFIND_${_l})
  ENDFOREACH(_l ${_libs})

      
  IF(_LFIND)

    MESSAGE( STATUS  "Trying to link mixed C/Fortran stuff...")

    SET(_CMAKELISTS ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB/CMakeLists.txt)

    # Create a CMakeLists.txt file which will generate the "flib" library

    # I
    FILE(WRITE ${_CMAKELISTS}
      "PROJECT(FortranTest C Fortran)\n"
      "SET(CMAKE_Fortran_FLAGS \"${TMP_Fortran_FLAGS}\")\n")

    # II
    FOREACH(_l ${_libs})
      IF(_LFIND_${_l})
        GET_FILENAME_COMPONENT(_LFINDDIR_l ${_LFIND_${_l}} PATH)
        FILE(APPEND ${_CMAKELISTS}
          "LINK_DIRECTORIES(${_LFINDDIR_l})\n")
      ENDIF(_LFIND_${_l})
    ENDFOREACH(_l ${_libs})

    # III
    FILE(APPEND ${_CMAKELISTS}
      "ADD_LIBRARY(flib ftest.f)\n"
      "ADD_EXECUTABLE(testF77libs testF77libs.c)\n")

    # IV Pfff...
    FOREACH(_l ${_libs})
      IF(_LFIND_${_l})
        FILE(APPEND ${_CMAKELISTS}
          "TARGET_LINK_LIBRARIES(testF77libs flib ${_l})\n")
      ENDIF(_LFIND_${_l})
    ENDFOREACH(_l ${_libs})


    # Use TRY_COMPILE to make the target "flib"
    TRY_COMPILE(FORT_LIBS_WORK
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB
      ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/FLIB
      flib
      OUTPUT_VARIABLE MY_OUTPUT)
    
    IF(FORT_LIBS_WORK)
      MESSAGE(STATUS "Ok, mixed C/Fortran can be linked with ${fc_libraries}")

      SET(FORTRAN_COMPILER_LIB_DIRECTORIES)
      SET(FORTRAN_LIBRARIES)
      FOREACH(_l ${_libs})
        IF(_LFIND_${_l})
          GET_FILENAME_COMPONENT(_LFINDDIR_l ${_LFIND_${_l}} PATH)
          LIST(APPEND FORTRAN_COMPILER_LIB_DIRECTORIES ${_LFINDDIR_l})
          LIST(APPEND FORTRAN_LIBRARIES ${_l})
        ENDIF(_LFIND_${_l})
      ENDFOREACH(_l ${_libs})
      
    ELSE(FORT_LIBS_WORK)
      MESSAGE(STATUS "No, mixed C/Fortran cannot be linked with ${fc_libraries}")
    ENDIF(FORT_LIBS_WORK)
  ENDIF(_LFIND)
    
  IF(FORT_LIBS_WORK)
    # the goal is to set this to something useful
    SET(FORTRAN_LIBRARIES ${FORTRAN_LIBRARIES} CACHE STRING "fortran libraries required to link with f90 generated code" FORCE)
    SET(FORTRAN_COMPILER_LIB_DIRECTORIES ${FORTRAN_COMPILER_LIB_DIRECTORIES} CACHE STRING "directories where fortran libraries can be found")
    MESSAGE(STATUS "Using FORTRAN_LIBRARIES ${FORTRAN_LIBRARIES}")
    MESSAGE(STATUS "Using FORTRAN_COMPILER_LIB_DIRECTORIES ${FORTRAN_COMPILER_LIB_DIRECTORIES}")
    
    SET(HAVE_FORTRAN_LIBRARIES TRUE CACHE INTERNAL "successfully found fortran libraries")
    SET(iopt ${imax})
    LINK_DIRECTORIES(${FORTRAN_COMPILER_LIB_DIRECTORIES})
  ELSE(FORT_LIBS_WORK)
    MATH(EXPR iopt ${iopt}+1)
  ENDIF(FORT_LIBS_WORK)
ENDWHILE(${iopt} LESS ${imax})
  
IF(NOT HAVE_FORTRAN_LIBRARIES)
  MESSAGE(FATAL_ERROR "Could not find valid fortran link libraries. Try setting FORTRAN_COMPILER_DIRECTORY or FORTRAN_LIBRARIES")
ENDIF(NOT HAVE_FORTRAN_LIBRARIES)

ENDIF(NOT HAVE_FORTRAN_LIBRARIES)
