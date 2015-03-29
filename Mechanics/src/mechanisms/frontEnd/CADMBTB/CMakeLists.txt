# Set minimum version for cmake
cmake_minimum_required(VERSION 2.8.7)

# MACOSX_RPATH required from cmake 3.0.
# See http://www.kitware.com/blog/home/post/510
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 OLD)
endif()

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../../CADMBTB ${CMAKE_CURRENT_SOURCE_DIR}/../TIMERMBTB)


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

  ELSE(OpenCASCADE_FOUND)
    MESSAGE(STATUS "OpenCascade Libraries not found in standard paths.")
  ENDIF(OpenCASCADE_FOUND)

endif(OCE_FOUND)

if(APPLE)
  find_package(PythonFull REQUIRED)
  include(FindPythonModule)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
else() # the previous instructions do not work on ubuntu 14.4 (and others)
       # the following instructions do not work on Yosemite.
IF(NOT WITH_PY3)
   # we still prefer python2 for now
   SET(PythonInterp_FIND_VERSION "2")
   SET(PythonInterp_FIND_VERSION_MAJOR "2")
   SET(PythonLibs_FIND_VERSION "2")
ELSE(NOT WITH_PY3)
  SET(PythonInterp_FIND_VERSION "3")
  SET(PythonInterp_FIND_VERSION_MAJOR "3")
  SET(PythonLibs_FIND_VERSION "3")
ENDIF(NOT WITH_PY3)

FIND_PACKAGE(PythonInterp)

# should world for all recent cmake version
# setting PythonInterp_FIND_VERSION to the version string is borken with cmake 3
SET(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

FIND_PACKAGE(PythonLibs)
endif()

MESSAGE("Python libs version : ${PYTHONLIBS_VERSION_STRING}")

MESSAGE("PYTHON EXECUTABLE : ${PYTHON_EXECUTABLE}")
MESSAGE("Python_ADDITIONAL_VERSIONS : ${Python_ADDITIONAL_VERSIONS}")
MESSAGE("Python_FRAMEWORKS : ${Python_FRAMEWORKS}")

MESSAGE("PYTHON INCLUDE DIR : ${PYTHON_INCLUDE_DIR}")
MESSAGE("PYTHON LIBRARIES : ${PYTHON_LIBRARIES}")




INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

# get python dist-packages
# warning python may not know cmake_install_prefix
# installation in standard python dist package with:
# -DSTD_PYTHON_DIST_PACKAGES=TRUE 
#
# installation in specific place:
# -DPYTHON_DIST_PACKAGES=/my/specific/place

IF(STD_PYTHON_DIST_PACKAGES)
  EXECUTE_PROCESS ( COMMAND 
    ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE PYTHON_DIST_PACKAGES 
    OUTPUT_STRIP_TRAILING_WHITESPACE)
ELSE(STD_PYTHON_DIST_PACKAGES)
  IF(NOT PYTHON_DIST_PACKAGES)
    EXECUTE_PROCESS ( COMMAND 
      ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib;print(get_python_lib(0,0,\"\"))"
      OUTPUT_VARIABLE PYTHON_DIST_PACKAGES 
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ENDIF(NOT PYTHON_DIST_PACKAGES)
ENDIF(STD_PYTHON_DIST_PACKAGES)

MESSAGE(STATUS "PYTHON_DIST_PACKAGES : ${PYTHON_DIST_PACKAGES}")

# do not skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
ASSERT(CMAKE_INSTALL_LIBDIR)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)



COMPILE_WITH(Boost REQUIRED)




SET(LIB_FRONTEND_CADMBTB_DEP CADMBTB TIMERMBTB)


INCLUDE(SiconosSwigHelper)

SICONOS_SWIG_MODULE(cadmbtb "cadmbtb;")
SWIG_LINK_LIBRARIES(cadmbtb ${PYTHON_LIBRARIES} ${LIB_FRONTEND_CADMBTB_DEP}
  ${OpenCASCADE_LIBRARIES})

