# It is the library libCADMBTB.so providing the CAD services.
# Require openCascade includes and libraries


SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/)

INCLUDE(SiconosTools)

APPEND_CXX_FLAGS("-DHAVE_CONFIG_H")

SET(OCE_TOOLKITS "TKernel"  "TKMath" "TKService" "TKV3d" "TKV2d" "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")

#find_package(OCE COMPONENTS ${OCE_TOOLKITS})

if(OCE_FOUND)
    message(STATUS "Found OCE version ${OCE_VERSION}")
    if(NOT OCE_ALL_FOUND)
      set(OCE_FOUND false)
      message(WARNING "Ignoring OCE installation due to missing toolkit(s): ${OCE_MISSING_TOOLKITS}")
    endif(NOT OCE_ALL_FOUND)
endif(OCE_FOUND)

if(OCE_FOUND)
  # Include files reside in ${OCE_INCLUDE_DIRS};
  include_directories(${OCE_INCLUDE_DIRS})
  LINK_LIBRARIES(CADMBTB ${OCE_LIBRARIES} TIMERMBTB)
  # We do not need library path, they will be automatically imported.
else(OCE_FOUND)
  # OCE not found; either it is not found and user
  # has to set OCE_DIR to the directory containing
  # OCEConfig.cmake, or OCE is not installed and we
  # try to find OpenCascade files.
  message(STATUS "OCE not found.  Try to find OpenCascade files.")
  
  FIND_PACKAGE(OpenCASCADE REQUIRED COMPONENTS ${OCE_TOOLKITS})

  IF(OpenCASCADE_FOUND) 
    MESSAGE(STATUS "OpenCascade Libraries found")
    message(STATUS "OpenCASCADE_INCLUDE_DIR = " ${OpenCASCADE_INCLUDE_DIR})
    message(STATUS "OpenCASCADE_LIBRARIES = " ${OpenCASCADE_LIBRARIES})
    message(STATUS "OpenCASCADE_LINK_DIRECTORY = " ${OpenCASCADE_LINK_DIRECTORY})
    include_directories(${OpenCASCADE_INCLUDE_DIR})
    LINK_LIBRARIES(CADMBTB ${OpenCASCADE_LIBRARIES}  TIMERMBTB)
  ELSE(OpenCASCADE_FOUND)
    MESSAGE(STATUS "OpenCascade Libraries not found in standard paths.")
  ENDIF(OpenCASCADE_FOUND)

endif(OCE_FOUND)

#AUX_SOURCE_DIRECTORY(src ${CMAKE_CURRENT_SOURCE_DIR}/pathOCC/)

include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../TIMERMBTB )
COMPILE_WITH(Boost REQUIRED)

enable_language (Fortran)

IF(APPLE)
SET (CMAKE_OSX_ARCHITECTURES "x86_64" )
ENDIF(APPLE)

ADD_LIBRARY(CADMBTB SHARED CADMBTB_API.cpp CADMBTB_PYTHON_API.cpp CADMBTB_internalTools.cpp myBRepExtrema_ExtFF.cxx myExtrema_ExtSS.cxx myExtrema_GenExtSS.cxx myExtrema_FuncExtSS.cxx mymath_FunctionSetRoot.cxx qnb.f) 
