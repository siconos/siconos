#
# Some convenience macros
#

# -- Basic list manipulation --
# Get first element of list var
MACRO(CAR var)
  SET(${var} ${ARGV1})
ENDMACRO(CAR)

# get elements in list var minus the first one.
MACRO(CDR var junk)
  SET(${var} ${ARGN})
ENDMACRO(CDR)

# LIST(APPEND ...) is not correct on <COMPILER>_FLAGS 
MACRO(APPEND_FLAGS)
  CAR(_V ${ARGV})
  CDR(_F ${ARGV})
  SET(${_V} "${${_V}} ${_F}")
ENDMACRO(APPEND_FLAGS)

# The use of ADD_DEFINITION results in a warning with Fortran compiler
MACRO(APPEND_C_FLAGS)
  APPEND_FLAGS(CMAKE_C_FLAGS ${ARGV})
ENDMACRO(APPEND_C_FLAGS)

MACRO(APPEND_CXX_FLAGS)
  APPEND_FLAGS(CMAKE_CXX_FLAGS ${ARGV})
ENDMACRO(APPEND_CXX_FLAGS)

MACRO(APPEND_Fortran_FLAGS)
  APPEND_FLAGS(CMAKE_Fortran_FLAGS ${ARGV})
ENDMACRO(APPEND_Fortran_FLAGS)

# Get a list of uninitialized cache variables and append
# them to vars.
macro(get_uninitialized_vars vars)
  get_cmake_property(CACHE_VARS CACHE_VARIABLES)
  foreach(CACHE_VAR ${CACHE_VARS})
    get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
    if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
      list(APPEND ${vars} -D${CACHE_VAR}=${${CACHE_VAR}})
    endif()
  endforeach()
endmacro()

# Scans DIRS (list of directories) and returns a list of all files in those dirs
# matching extensions defined in SRC_EXTS list.
# Results are saved in SOURCES_FILES
#
# Usage:
# set(src_dirs dir1 dir2)
# get_sources(src_dirs)
macro(get_sources)
  set(SOURCES_FILES)
  foreach(DIR ${ARGV})
    foreach(_EXT ${SRC_EXTS})
      file(GLOB FILES_LIST ${DIR}/*.${_EXT})
      if(FILES_LIST)
	list(APPEND SOURCES_FILES ${FILES_LIST})
      endif()
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES SOURCES_FILES)
endmacro()

# Scans DIRS (list of directories) and returns a list of all files in those dirs
# matching extensions defined in HDR_EXTS list.
# Results are saved in HDRS_FILES
#
# Usage:
# set(src_dirs dir1 dir2)
# get_headers(src_dirs)
macro(get_headers DIRS)
  set(HDRS_FILES)
  foreach(DIR ${ARGV})
    foreach(_EXT ${HDR_EXTS})
      file(GLOB FILES_LIST ${DIR}/*.${_EXT})
      if(FILES_LIST)
	list(APPEND HDRS_FILES ${FILES_LIST})
      endif()
    endforeach()
  endforeach()
  list(REMOVE_DUPLICATES HDRS_FILES)
endmacro()

# -- returns a list of source files extension --
# Results in var ALL_EXTS
macro(get_standard_ext)
  set(ALL_EXTS)
  foreach(_EXT
      ${CMAKE_CXX_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_C_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_Java_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_RC_SOURCE_FILE_EXTENSIONS})
    list(APPEND ALL_EXTS ${_EXT})
  endforeach()
  list(REMOVE_DUPLICATES ALL_EXTS)
endmacro()

# Print cmake variable 'V' value
MACRO(PRINT_VAR V)
  MESSAGE(STATUS "${V} = ${${V}}")
ENDMACRO(PRINT_VAR V)


# =======================================
# For a given package name, try to find
# corresponding headers and libraries and
# add them to the include directories
# and list of linked libraries.
#
# It sets (if found):
# - SICONOS_INCLUDE_DIRECTORIES with the list
# of directories of headers required for siconos to work with
# - SICONOS_LINK_LIBRARIES with the list of external libraries
# (full path!) needed by siconos project.
#
# Usage :
#  compile_with(Packagename options)
#
# with the same 'options' as find_package
# (see http://www.cmake.org/cmake/help/v3.0/command/find_package.html?highlight=find_package)
MACRO(COMPILE_WITH)
  # Get package name and extra args ...
  CAR(_NAME ${ARGV})
  CDR(_REST ${ARGV})
  CAR(_REQ ${_REST})
  IF(_REST)
    CDR(_RREST ${_REST})
    CAR(_COMP ${_RREST})
  ENDIF(_REST)
  STRING(TOUPPER ${_NAME} _UNAME)
  LIST(APPEND _NAMES ${_NAME})
  LIST(APPEND _NAMES ${_UNAME})
  SET(_FOUND)
  FIND_PACKAGE(${ARGV})
  
  FOREACH(_N ${_NAMES})
    IF(${_N}_FOUND)
      SET(_FOUND TRUE)
      # add headers dirs into 'include' path
      # INCLUDE_DIR var name depends on FindNAME
      # We try to check the standard var names.
      if(DEFINED ${_N}_INCLUDE_DIRS)
	remember_include_directories("${${_N}_INCLUDE_DIRS}")
      endif()
      if(DEFINED ${_N}_INCLUDE_DIR)
	remember_include_directories("${${_N}_INCLUDE_DIR}")
      endif()
      if(DEFINED ${_N}_INCLUDE_PATH)
	remember_include_directories("${${_N}_INCLUDE_DIR}")
      endif()
      # Now we set list of libs that must be linked with.
      if(DEFINED ${_N}_LIBRARIES)
	set(SICONOS_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES}
	  ${${_N}_LIBRARIES})
      endif()
      # And the compiler flags
      if(DEFINED ${_N}_DEFINITIONS)
        APPEND_C_FLAGS(${${_N}_DEFINITIONS})
	APPEND_CXX_FLAGS(${${_N}_DEFINITIONS})
      endif()
    endif()
  endforeach()
  list(REMOVE_DUPLICATES SICONOS_LINK_LIBRARIES)
  set(SICONOS_LINK_LIBRARIES ${SICONOS_LINK_LIBRARIES}
    ${${_N}_LIBRARIES} CACHE INTERNAL "List of external libraries.")
  set(_N)
  set(_NAME) 
  set(_UNAME)
  set(_NAMES)
ENDMACRO(COMPILE_WITH)

# ==== Save directories required for include_directory ===
# 
# Set variable SICONOS_INCLUDE_DIRECTORIES with the list
# of directories of headers required for siconos to work with
# its dependencies.
# Usage :
# set(dirs d1 d2 d3)
# remember_include_directories(${dirs})
#  --> save d1, d2, d3 into SICONOS_INCLUDE_DIRECTORIES
# 
MACRO(REMEMBER_INCLUDE_DIRECTORIES _DIRS)
  FOREACH(_D ${_DIRS})
    LIST(APPEND SICONOS_INCLUDE_DIRECTORIES ${_D})
  ENDFOREACH()
  list(REMOVE_DUPLICATES SICONOS_INCLUDE_DIRECTORIES)
  set(SICONOS_INCLUDE_DIRECTORIES ${SICONOS_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "Include directories for external dependencies.")

ENDMACRO()

# ==== Save directories required for include_directory ===
# 
# Set variable ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES with the list
# of directories of headers of each siconos component.
#
# Usage :
# set(dirs d1 d2 d3)
# remember_local_include(${dirs})
#  --> save d1, d2, d3 into ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
#
# mind the ${CMAKE_CURRENT_SOURCE_DIR} below!
macro(remember_local_include_directories _DIRS)
  foreach(_D ${_DIRS})
    list(APPEND ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
      ${CMAKE_CURRENT_SOURCE_DIR}/${_D})
  endforeach()
  list(REMOVE_DUPLICATES ${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES)
  set(${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES
    ${${PROJECT_NAME}_LOCAL_INCLUDE_DIRECTORIES}
    CACHE INTERNAL "Include directories for external dependencies.")
endmacro()


MACRO(WRITE_NOTES)
  IF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
    # a note file for the dashboard
    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Notes)
    FILE(WRITE ${CMAKE_BINARY_DIR}/Testing/Notes/Build "git sha1 : ${SOURCE_ABBREV_GIT_SHA1}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "cmake version : ${CMAKE_VERSION}\n")
    # the default buildname
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "System name : ${CMAKE_SYSTEM_NAME}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Processor   : ${CMAKE_SYSTEM_PROCESSOR}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler : ${CMAKE_C_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler version : ${C_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler : ${CMAKE_CXX_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler version : ${CXX_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler : ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler version : ${Fortran_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "BLAS libraries : ${BLAS_LIBRARIES}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "LAPACK libraries : ${LAPACK_LIBRARIES}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "all libraries : ${SICONOS_LINK_LIBRARIES}\n")
  ENDIF(IS_DIRECTORY ${CMAKE_BINARY_DIR}/Testing)
ENDMACRO(WRITE_NOTES)

MACRO(ASSERT VAR)
  IF (NOT DEFINED ${VAR})
    MESSAGE( FATAL_ERROR "ASSERTION ERROR : ${VAR} UNSET" )
  ENDIF()
ENDMACRO()    


# -------------------------------
# Set WITH_COMPONENT_OPT value
# depending on WITH_OPT value
# and the -D... entries.
#
# Example :
# cmake -DWITH_DOCUMENTATION = ON
# will activate all WITH_component_DOCUMENTATION for enabled components.
# while
# cmake -DWITH_kernel_DOCUMENTATION=ON
# will set WITH_DOCUMENTATION=ON and WITH_other_components=OFF
# 
# This will work (I hope ...) in standard cases but will probably
# failed after several cmake . with schizophrenic options
# like
# cmake -DWITH_kernel_DOCUMENTATION=ON path_to_srcs
# cmake -DWITH_DOCUMENTATION=OFF .
# In that case, user needs to reset all WITH_component_OPT.
#
# -------------------------------
macro(init_to_default_option OPT)
  # Each "WITH_component_OPT" is set to default value == WITH_OPT value.
  foreach(comp ${COMPONENTS_NAMES})
    if(NOT WITH_${comp}_${OPT})
      set(WITH_${comp}_${OPT} ${WITH_${OPT}} CACHE BOOL "initialize ${OPT} for component ${comp}.")
    endif()
    # We don't want to see all with_comp_opt in the GUI.
    mark_as_advanced(WITH_${comp}_${OPT})
 endforeach()

 # If one with_comp_opt is on, global with_opt must also be on
 foreach(comp ${COMPONENTS_NAMES})
   if(WITH_${comp}_${OPT})
     set(WITH_${OPT} ON  CACHE BOOL "initialize ${OPT}." FORCE)
     break()
   endif()
 endforeach()
endmacro()
