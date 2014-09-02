# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
cmake_policy(SET CMP0042 OLD)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../../CADMBTB ${CMAKE_CURRENT_SOURCE_DIR}/../TIMERMBTB)
SET_SOURCE_FILES_PROPERTIES(cadmbtb.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(cadmbtb python cadmbtb.i)


message(STATUS "Searching for OCE ....")
SET(OCE_TOOLKITS "TKernel" "PTKernel" "TKMath" "TKService" "TKV3d"  "TKBRep" "TKIGES" "TKSTL" "TKVRML" "TKSTEP" "TKSTEPAttr" "TKSTEP209" "TKSTEPBase" "TKShapeSchema" "TKGeomBase" "TKGeomAlgo" "TKG3d" "TKG2d" "TKXSBase" "TKPShape" "TKShHealing" "TKHLR" "TKTopAlgo" "TKMesh" "TKPrim" "TKCDF" "TKBool" "TKBO" "TKFillet" "TKOffset")

find_package(OCE COMPONENTS ${OCE_TOOLKITS})

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
  SWIG_LINK_LIBRARIES(cadmbtb ${OpenCASCADE_LIBRARIES} )
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

    SWIG_LINK_LIBRARIES(cadmbtb ${OpenCASCADE_LIBRARIES} )
  ELSE(OpenCASCADE_FOUND)
    MESSAGE(STATUS "OpenCascade Libraries not found in standard paths.")
  ENDIF(OpenCASCADE_FOUND)

endif(OCE_FOUND)

SET(Python_ADDITIONAL_VERSIONS "2.7;2.6")
FIND_PACKAGE(PythonInterp)

# on ubuntu precise with cmake >= 2.8.8 pythonlibs find python 3.0
# before python 2.7 (in Release but not in Debug ...)
#SET(Python_ADDITIONAL_VERSIONS "2.6;2.7")

FIND_PACKAGE(PythonLibs)

MESSAGE("Python libs version : ${PYTHONLIBS_VERSION_STRING}")

MESSAGE("PYTHON EXECUTABLE : ${PYTHON_EXECUTABLE}")
MESSAGE("Python_ADDITIONAL_VERSIONS : ${Python_ADDITIONAL_VERSIONS}")
MESSAGE("Python_FRAMEWORKS : ${Python_FRAMEWORKS}")

MESSAGE("PYTHON INCLUDE DIR : ${PYTHON_INCLUDE_DIR}")
MESSAGE("PYTHON LIBRARIES : ${PYTHON_LIBRARIES}")



INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})


COMPILE_WITH(Boost REQUIRED)




SET(LIB_FRONTEND_CADMBTB_DEP CADMBTB TIMERMBTB)

SWIG_LINK_LIBRARIES(cadmbtb ${PYTHON_LIBRARIES} ${LIB_FRONTEND_CADMBTB_DEP})


