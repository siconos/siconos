#
# Out of sources made mandatory
#
OPTION(IN_SOURCE_BUILD "if you really want a build in the source directory" OFF)

IF(NOT IN_SOURCE_BUILD)
  IF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    MESSAGE(FATAL_ERROR "CMake generation for this library is not allowed within the source directory! 
Remove the CMakeCache.txt file and try again from another folder, e.g.: 

   rm CMakeCache.txt 
   cd <somewhere (preferably a local place on your computer and not a network folder)>
   cmake <source directory of Numerics>

If you really need an in source build, then run : cmake -DIN_SOURCE_BUILD=ON
")
  ENDIF(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
ENDIF(NOT IN_SOURCE_BUILD)
