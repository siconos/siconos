SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})


FIND_PACKAGE(SiconosNumerics)
FIND_PACKAGE(SiconosKernel)

# #FIND_PACKAGE(OpenCascade)

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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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
COMPILE_WITH(Python_Numpy REQUIRED)



SET(MBTB_LINK_LIBRARIES TIMERMBTB CADMBTB MBTB ${SiconosKernel_LIBRARIES} ${SiconosNumerics_LIBRARIES} ${PYTHON_LIBRARIES})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MBTB ${SiconosKernel_INCLUDE_DIRS}
  ${SiconosNumerics_INCLUDE_DIRS} ${ATLAS_INCLUDE_DIR})

INCLUDE(SiconosSwigHelper)

SICONOS_SWIG_MODULE(mbtb "mbtb;")
target_link_libraries(${SWIG_MODULE_mbtb_REAL_NAME} ${MBTB_LINK_LIBRARIES})

