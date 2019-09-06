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

# Collect source files.
# 
# Usage:
#
# get_sources(<COMPONENT> DIRS <dirs list> EXCLUDE <files list>)
#
# Result : set (parent scope) <COMPONENT>_SRCS with files in
# dir1, dir2 matching standard extension for C,C++ and Fortran.
# Do not include files listed after EXCLUDE option.
#
# Remarks:
# - dir1, dir2 ... are relative to CMAKE_CURRENT_SOURCE_DIR
# 
function(get_sources COMPONENT)
  
  set(multiValueArgs DIRS EXCLUDE)
  cmake_parse_arguments(source "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Get list of extensions to be taken into account
  foreach(_EXT
      ${CMAKE_CXX_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_C_SOURCE_FILE_EXTENSIONS}
      ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS})
    list(APPEND SRC_EXTS ${_EXT})
  endforeach()
  list(REMOVE_DUPLICATES SRC_EXTS)

  # Scan all dirs and check all exts ...
  foreach(DIR IN LISTS source_DIRS)
    foreach(_EXT IN LISTS SRC_EXTS)
      if(${CMAKE_VERSION} VERSION_GREATER "3.12.0")
        file(GLOB FILES_LIST
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CONFIGURE_DEPENDS
          ${DIR}/*.${_EXT})
      else()
        file(GLOB FILES_LIST RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${DIR}/*.${_EXT})
      endif()
      if(FILES_LIST)
	list(APPEND SOURCES_FILES ${FILES_LIST})
      endif()
    endforeach()
  endforeach()
  if(SOURCES_FILES)
    list(LENGTH SOURCES_FILES _SOURCES_FILES_LEN)
    if (_SOURCES_FILES_LEN GREATER 1)
      list(REMOVE_DUPLICATES SOURCES_FILES)
    endif()
  endif()

  # Check if some sources are to be excluded from build
  foreach(_FILE IN LISTS source_EXCLUDE)
    if(${CMAKE_VERSION} VERSION_GREATER "3.12.0")
      file(GLOB _GFILE CONFIGURE_DEPENDS ${_FILE})
    else()
      file(GLOB _GFILE ${_FILE})
    endif()
    
    if(_GFILE)
      list(REMOVE_ITEM SOURCES_FILES ${_GFILE})
    else()
      message(WARNING "file to be excluded NOT FOUND : ${_FILE}")
    endif()
  endforeach()
  
  set(${COMPONENT}_SRCS ${SOURCES_FILES} PARENT_SCOPE)
endfunction()

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
      if (INSTALL_INTERNAL_HEADERS AND FILES_LIST)
	list(APPEND HDRS_FILES ${FILES_LIST})
      else()
        # filter out header paths containing the word "internal"
        # (stemming from component dir.. otherwise we'd have trouble
        # if workdir path contains the string "internal")
        foreach(_HDR ${FILES_LIST})
          if (_HDR AND NOT "${_HDR}" MATCHES "${_COMPONENT}/.*internal")
	    list(APPEND HDRS_FILES ${_HDR})
          endif()
        endforeach()
      endif()
    endforeach()
  endforeach()
  list(LENGTH HDRS_FILES _HDRS_FILES_LEN)
  if (_HDRS_FILES_LEN GREATER 1)
    list(REMOVE_DUPLICATES HDRS_FILES)
  endif()
endmacro()


# Print cmake variable 'V' value
macro(PRINT_VAR V)
  message(STATUS "${V} = ${${V}}")
endmacro()


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
    site_name(_SITE_NAME)
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Site Name: ${_SITE_NAME}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Processor   : ${CMAKE_SYSTEM_PROCESSOR}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler : ${CMAKE_C_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "C compiler version : ${CMAKE_C_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler : ${CMAKE_CXX_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "CXX compiler version : ${CMAKE_CXX_COMPILER_VERSION}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler : ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_BINARY_DIR}/Testing/Notes/Build "Fortran compiler version : ${CMAKE_Fortran_COMPILER_VERSION}\n")
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
  foreach(comp ${COMPONENTS})
    if(NOT WITH_${comp}_${OPT})
      set(WITH_${comp}_${OPT} ${WITH_${OPT}} CACHE BOOL "initialize ${OPT} for component ${comp}.")
    endif()
    # We don't want to see all with_comp_opt in the GUI.
    mark_as_advanced(WITH_${comp}_${OPT})
 endforeach()

 # If one with_comp_opt is on, global with_opt must also be on
 foreach(comp ${COMPONENTS})
   if(WITH_${comp}_${OPT})
     set(WITH_${OPT} ON  CACHE BOOL "initialize ${OPT}." FORCE)
     break()
   endif()
 endforeach()
endmacro()

# Display MPI search results
function(print_mpi_info lang)
  message("\n--------------------------- MPI ${lang} config ---------------------------")
  message("- compiler: ${MPI_${lang}_COMPILER}")
  message("- compile flags: ${MPI_${lang}_COMPILE_FLAGS}")
  message("- include path: ${MPI_${lang}_INCLUDE_PATH}")
  message("- link flags: ${MPI_${lang}_LINK_FLAGS}")
  message("- libraries: ${MPI_${lang}_LIBRARIES}")
  message("-------------------------------------------------------------------------\n")
endfunction()
