add_subdirectory(TIMERMBTB)
add_subdirectory(CADMBTB)
add_subdirectory(MBTB)
add_subdirectory(plugin/C60)
add_subdirectory(plugin/SliderCrank)
add_subdirectory(frontEnd)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

INCLUDE(SiconosTools)
INCLUDE(UseDoxygen OPTIONAL)

# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 

FIND_PACKAGE(Boost COMPONENTS graph)

IF(NOT Boost_GRAPH_FOUND)
  MESSAGE(STATUS "WARNING boost graph library not found.")
ENDIF(NOT Boost_GRAPH_FOUND)


enable_language (Fortran)

IF(APPLE)
SET (CMAKE_OSX_ARCHITECTURES "x86_64" )
ENDIF(APPLE)
