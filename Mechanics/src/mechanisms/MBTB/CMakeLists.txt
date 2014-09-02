# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
cmake_policy(SET CMP0042 OLD)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/)
FIND_PACKAGE(SiconosNumerics REQUIRED)
FIND_PACKAGE(SiconosKernel REQUIRED)
 

INCLUDE_DIRECTORIES(../../joints)
COMPILE_WITH(BLAS REQUIRED)

COMPILE_WITH(LAPACK REQUIRED COMPLETE)




##FIND_PACKAGE(ATLAS REQUIRED)

INCLUDE(UseDoxygen OPTIONAL)

include_directories ( ${SiconosMechanics_INCLUDE_DIRS} ${SiconosKernel_INCLUDE_DIRS} ${SiconosNumerics_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../CADMBTB ${CMAKE_CURRENT_SOURCE_DIR}/../TIMERMBTB $ ${CBLAS_INCLUDE_PATH} )

COMPILE_WITH(Boost REQUIRED)
# FIND_PACKAGE( ... COMPONENTS <component>) does not mention <component> in failure message!
# a workaround 

FIND_PACKAGE(Boost COMPONENTS graph)

IF(NOT Boost_GRAPH_FOUND)
  MESSAGE(STATUS "WARNING boost graph library not found.")
ENDIF(NOT Boost_GRAPH_FOUND)


#LINK_LIBRARIES(MBTB ${CMAKE_CURRENT_BINARY_DIR}/../CADMBTB/libCADMBTB.so ${CMAKE_CURRENT_BINARY_DIR}/../TIMERMBTB/libTIMERMBTB.so)

LINK_LIBRARIES(MBTB TIMERMBTB CADMBTB ${SiconosMechanics_LIBRARY} ${SiconosKernel_LIBRARIES} ${SiconosNumerics_LIBRARIES})

add_library (MBTB SHARED MBTB_PYTHON_API.cpp MBTB_Body.cpp MBTB_internalTool.cpp MBTB_DATA.cpp MBTB_Contact.cpp MBTB_ContactRelation.cpp MBTB_FC3DContactRelation.cpp MBTB_TimeStepping.cpp MBTB_TimeSteppingProj.cpp MBTB_TimeSteppingCombinedProj.cpp MBTB_MoreauJeanOSI.cpp MBTB_JointR.cpp)
