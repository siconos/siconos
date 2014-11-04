# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()

FIND_PACKAGE(Boost COMPONENTS graph)

FIND_PACKAGE(SiconosMechanics REQUIRED)

IF(NOT Boost_GRAPH_FOUND)
  MESSAGE(STATUS "WARNING boost graph library not found.")
ENDIF(NOT Boost_GRAPH_FOUND)



add_subdirectory(TIMERMBTB)
add_subdirectory(CADMBTB)
add_subdirectory(MBTB)
add_subdirectory(plugin/C60)
add_subdirectory(plugin/SliderCrank)
add_subdirectory(frontEnd)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
INCLUDE(SiconosTools)
INCLUDE(UseDoxygen OPTIONAL)


# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 


# IF(APPLE)
# SET (CMAKE_OSX_ARCHITECTURES "x86_64" )
# ENDIF(APPLE)
