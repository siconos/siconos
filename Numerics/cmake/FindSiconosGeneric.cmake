# Find the Siconos ${_project} includes and libraries.
# The following variables are set if ${_project} is found.  If ${_project} is not
# found, Siconos${_project}_FOUND is set to false.
#  Siconos${_project}_FOUND        - True when the Siconos${_project} include directory is found.
#  Siconos${_project}_INCLUDE_DIRS - the path to where the Siconos ${_project} include files are.
#  Siconos${_project}_LIBRARY_DIRS - The path to where the Siconos library files are.
#  Siconos${_project}_LIBRARIES    - The libraries to link against Siconos ${_project}

# One may want to use a specific ${_project} Library by setting
# Siconos${_project}_LIBRARY_DIRECTORY before FIND_PACKAGE(Siconos${_project})
INCLUDE(FindPackageHandleStandardArgs)

MACRO(FIND_SICONOS_COMPONENT _project _header_file)

  # We should try to use Windows-GNU to some extend ... --xhub
  IF(CMAKE_SYSTEM_NAME MATCHES Windows)
    if (NOT MINGW)
      SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "" ${CMAKE_FIND_LIBRARY_PREFIXES})
      SET(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".dll.a" ".a")
    endif()
    SET(EXE_EXT ".exe")
  ELSE()
    SET(EXE_EXT)
  ENDIF()

  IF(Siconos${_project}_LIBRARY_DIRECTORY)
    MESSAGE(STATUS "Looking for ${_project} library in ${Siconos${_project}_LIBRARY_DIRECTORY}")
    FIND_LIBRARY(Siconos${_project}_LIBRARY Siconos${_project} PATHS "${Siconos${_project}_LIBRARY_DIRECTORY}" NO_DEFAULT_PATH)
    IF(Siconos${_project}_LIBRARY)
      MESSAGE(STATUS "Found : ${Siconos${_project}_LIBRARY}")
    ENDIF(Siconos${_project}_LIBRARY)
  ELSE(Siconos${_project}_LIBRARY_DIRECTORY)
    FIND_LIBRARY(Siconos${_project}_LIBRARY Siconos${_project} PATHS ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH)
  ENDIF(Siconos${_project}_LIBRARY_DIRECTORY)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(Siconos${_project} 
    REQUIRED_VARS Siconos${_project}_LIBRARY)

  # We should use the second version of FIND_PACKAGE_HANDLE_STANDARD_ARGS ... -- xhub
  string(TOUPPER ${_project} _PROJECT)
  SET(Siconos${_project}_FOUND ${SICONOS${_PROJECT}_FOUND})

  IF(Siconos${_project}_LIBRARY)
    SET(Siconos${_project}_LIBRARIES ${Siconos${_project}_LIBRARY})
    GET_FILENAME_COMPONENT(Siconos${_project}_LIBRARY_NAME ${Siconos${_project}_LIBRARY} NAME)
    GET_FILENAME_COMPONENT(Siconos${_project}_LIBRARY_DIRS ${Siconos${_project}_LIBRARY} PATH)
    GET_FILENAME_COMPONENT(Siconos${_project}_LIBRARY_DIRS_DIR ${Siconos${_project}_LIBRARY_DIRS} PATH)
    GET_FILENAME_COMPONENT(Siconos${_project}_LIBRARY_DIRS_DIR_DIR ${Siconos${_project}_LIBRARY_DIRS_DIR} PATH)

    # It looks like cmake does stupid things with the order in HINTS ...
    SET(Siconos${_project}_LOCATION_DIRS ${Siconos${_project}_LIBRARY_DIRS_DIR} ${Siconos${_project}_LIBRARY_DIRS_DIR_DIR})

    FOREACH(_DIR ${Siconos${_project}_LOCATION_DIRS})
      FIND_PATH(Siconos${_project}_INCLUDE_DIRS ${_header_file}
        HINTS ${_DIR}
        ENV PATH
        PATH_SUFFIXES include/Siconos/${_project})
    ENDFOREACH()

    IF(NOT Siconos${_project}_INCLUDE_DIRS)
      IF(Siconos${_project}_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR
          "Required Siconos ${_project} headers not found. Please specify headers location in CMAKE_INCLUDE_PATH")
      ENDIF(Siconos${_project}_FIND_REQUIRED)
    ELSE(NOT Siconos${_project}_INCLUDE_DIRS)
      MESSAGE(STATUS "Include directory for Siconos${_project}: ${Siconos${_project}_INCLUDE_DIRS}")
    ENDIF(NOT Siconos${_project}_INCLUDE_DIRS)

  ELSE(Siconos${_project}_LIBRARY)
    IF(Siconos${_project}_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR
        "Required Siconos ${_project} library not found. Please specify library location in Siconos${_project}_LIBRARY_DIRECTORY")
    ENDIF(Siconos${_project}_FIND_REQUIRED)
  ENDIF(Siconos${_project}_LIBRARY)

ENDMACRO()
