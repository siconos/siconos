# ===================================================
# Set main parameters for the Siconos cmake project
#
# ===================================================
include(SiconosTools)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif()

# =========== Windows stuff ... ===========
include(WindowsSiconosSetup)

# --------- CMake project internal variables ---------

# Siconos current version
include(SiconosVersion)

# File used to print tests setup messages.
set(TESTS_LOGFILE ${CMAKE_BINARY_DIR}/tests.log)
# get system architecture 
# https://raw.github.com/petroules/solar-cmake/master/TargetArch.cmake
include(TargetArch)
target_architecture(SYSTEM_ARCHITECTURE)
if(WITH_SYSTEM_INFO)
  include(CMakePrintSystemInformation)
  message(STATUS "SYSTEM ARCHITECTURE: ${SYSTEM_ARCHITECTURE}")
endif()

# extensions of headers files that must be taken into account
set(HDR_EXTS h hpp)

# To include or not unstable source files
set(WITH_UNSTABLE FALSE)

# dirs of 'local' headers. Must be filled by each component.
set(${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
  CACHE INTERNAL "Local include directories.")

set(SICONOS_INCLUDE_DIRECTORIES
  CACHE INTERNAL "Include directories for external dependencies.")

set(SICONOS_LINK_LIBRARIES ${CMAKE_DL_LIBS}
  CACHE INTERNAL "List of external libraries.")
# CMAKE_DL needed on Debian (and others?) for dlopen or equivalent (used for plugins)

set(${PROJECT_NAME}_LOCAL_LIBRARIES
  CACHE INTERNAL "List of siconos components libraries.")

set(installed_targets ${installed_targets}
  CACHE INTERNAL "List of installed libraries for the siconos project.")

if(CMAKE_VERSION VERSION_LESS 3)
  set(PRIVATE LINK_PRIVATE CACHE INTERNAL "")
  set(PUBLIC LINK_PUBLIC CACHE INTERNAL "")
else()
  set(PRIVATE PRIVATE CACHE INTERNAL "")
  set(PUBLIC PUBLIC CACHE INTERNAL "")
endif()

set(tests_timeout 15 CACHE INTERNAL "Limit time for tests (in seconds)")

# extensions of source files that must be taken into account
get_standard_ext()
set(SRC_EXTS ${ALL_EXTS})

if(WITH_GIT)
  find_package(Git)
  MESSAGE(STATUS "git executable : ${GIT_EXECUTABLE}")
  MESSAGE(STATUS "git command : ${GITCOMMAND}")
  MESSAGE(STATUS "git update options : ${GIT_UPDATE_OPTIONS}")
     
  SET(CTEST_GIT_COMMAND "${GIT_EXECUTABLE}" )
  SET(UPDATE_COMMAND "${GITCOMMAND}")
  SET(UPDATE_OPTIONS "${GIT_UPDATE_OPTIONS}")
      
  EXECUTE_PROCESS(COMMAND 
    ${GIT_EXECUTABLE} log -n 1 --pretty=format:%h 
    OUTPUT_VARIABLE SOURCE_ABBREV_GIT_SHA1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  
  EXECUTE_PROCESS(COMMAND 
    ${GIT_EXECUTABLE} log -n 1 --pretty=format:%H
    OUTPUT_VARIABLE SOURCE_GIT_SHA1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# ---- Python ---
# (interp and lib)
if(WITH_PYTHON_WRAPPER OR WITH_DOCUMENTATION)
  find_package(PythonFull REQUIRED)
  include(FindPythonModule)
  get_filename_component(PYTHON_EXE_NAME ${PYTHON_EXECUTABLE} NAME)

  # --- xml schema. Used in tests. ---
  if(WITH_XML)
    set(SICONOS_XML_SCHEMA "${CMAKE_SOURCE_DIR}/config/xmlschema/SiconosModelSchema-V3.7.xsd")
    IF(NOT NO_RUNTIME_BUILD_DEP)
      find_python_module(lxml REQUIRED)
    ENDIF(NOT NO_RUNTIME_BUILD_DEP)
  endif()

endif()

IF(NOT SIZE_OF_CSI)
  INCLUDE(CheckTypeSize)
  CHECK_TYPE_SIZE("size_t" SIZE_OF_CSI)
ENDIF(NOT SIZE_OF_CSI)

IF ("${SIZE_OF_CSI}" EQUAL 8)
  SET(SICONOS_INT64 TRUE)
ELSE ("${SIZE_OF_CSI}" EQUAL 8)
  SET(SICONOS_INT64 FALSE)
ENDIF ("${SIZE_OF_CSI}" EQUAL 8)

# =========== install setup ===========

# Install lib directory 32, 64 etc. on Fedora, Debian 
# http://public.kitware.com/Bug/view.php?id=11964
# See also http://www.cmake.org/cmake/help/v3.0/module/GNUInstallDirs.html?highlight=gnuinstalldirs
include(GNUInstallDirs)
# Set prefix path for libraries installation
# --> means that any library target will be installed
# in CMAKE_INSTALL_PREFIX/_install_lib
if(${PROJECT_NAME}_INSTALL_LIB_DIR)
  set(_install_lib ${${PROJECT_NAME}_INSTALL_LIB_DIR})
else()
  ASSERT(CMAKE_INSTALL_LIBDIR)
  set(_install_lib ${CMAKE_INSTALL_LIBDIR})
  set(${PROJECT_NAME}_INSTALL_LIB_DIR ${_install_lib})
endif()

# install cmake macros
# ==> those that may be required during call to siconos script.
set(cmake_macros
  SiconosTools.cmake
  FindQGLViewer.cmake
  )
foreach(file ${cmake_macros})
  install(FILES cmake/${file} DESTINATION share/${PROJECT_NAME}/cmake)
endforeach()

# =========== uninstall target ===========
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  echo >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
  COMMAND cat ${CMAKE_CURRENT_BINARY_DIR}/python_install_manifest.txt >> ${CMAKE_CURRENT_BINARY_DIR}/install_manifest.txt
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# # =========== RPATH stuff ===========
# # we follow recommendation of https://cmake.org/Wiki/CMake_RPATH_handling

if(APPLE)
  set(CMAKE_MACOSX_RPATH TRUE)

# The following settings were copied from
# https://cmake.org/Wiki/CMake_RPATH_handling
# to avoid the rpath issue that appears on OS X El Capitan

  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) # Changed to TRUE by A.O.

  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")
else(APPLE)
  # =========== RPATH stuff ===========
  # do not skip the full RPATH for the build tree
  if(FORCE_SKIP_RPATH)
    set(CMAKE_SKIP_BUILD_RPATH TRUE)
  else(FORCE_SKIP_RPATH)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
  endif(FORCE_SKIP_RPATH)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

  # when building a binary package, it makes no sense to add this rapth
  if(NOT FORCE_SKIP_RPATH)
    # the RPATH to be used when installing
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  endif(NOT FORCE_SKIP_RPATH)

  # don't add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif(APPLE)

# init all common options for enabled components
set(common_options DOCUMENTATION TESTING UNSTABLE PYTHON_WRAPPER
  DOXYGEN_WARNINGS DOXY2SWIG)
foreach(opt ${common_options})
  init_to_default_option(${opt})
endforeach()
# ========= Documentation =========
if(WITH_DOCUMENTATION OR WITH_DOXY2SWIG OR WITH_DOXYGEN_WARNINGS)
  set(USE_DOXYGEN TRUE)
endif()

if(WITH_DOCUMENTATION)
  # temporary option?
  set(USE_SPHINX TRUE)
endif()

