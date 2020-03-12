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

# Print cmake variable 'V' value
macro(PRINT_VAR V)
  message(STATUS "${V} = ${${V}}")
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


# Try to provide some hints for a find_package call.
# 
# Usage:
#
# set_find_package_hints(NAME <name> MODULE <mod>)
#
# 
# Result : set (parent scope) _<NAME>_SEARCH_OPTS and _<NAME>_INC_SEARCH_OPTS
# that can be used in find_path (INC_SEARCH) and find_library calls.
#
# These variables are filled either with <NAME>_ROOT value if set
# or using pkg-config information, if available.
#
# See examples of use in FindCPPUNIT.cmake or FindSuperLU.cmake.
# 
function(set_find_package_hints)
  set(oneValueArgs NAME MODULE)

  cmake_parse_arguments(pkg "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(${pkg_NAME}_ROOT)
    set(_${pkg_NAME}_SEARCH_OPTS
      HINTS ${${pkg_NAME}_ROOT} NO_DEFAULT_PATH PARENT_SCOPE)
    set(_${pkg_NAME}_INC_SEARCH_OPTS
      HINTS ${${pkg_NAME}_ROOT} NO_DEFAULT_PATH PARENT_SCOPE)
  else()
    # Try pkgconfig
    find_package(PkgConfig QUIET)
    pkg_check_modules(PKGC_${pkg_NAME} ${pkg_MODULE} QUIET)
    if(PKGC_${pkg_NAME}_FOUND)
      set(_${pkg_NAME}_INC_SEARCH_OPTS "HINTS ${PKGC_${pkg_NAME}_INCLUDE_DIRS}"
        PARENT_SCOPE)
    endif()
    set(_${pkg_NAME}_SEARCH_OPTS
      HINTS ${PKGC_${pkg_NAME}_LIBRARY_DIRS} ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
      PARENT_SCOPE)
  endif()
  
endfunction()

# ------------------------------------
# Get the list of subdirectories
# of a given dir
# Useful only in examples ...
# ------------------------------------
macro(get_subdirectories result current_dir)
  file(GLOB subdirs RELATIVE ${current_dir} ${current_dir}/*)
  set(dirs "")
  foreach(_dir ${subdirs})
    if(IS_DIRECTORY ${current_dir}/${_dir})
      list(APPEND dirs ${_dir})
    endif()
  endforeach()
  set(${result} ${dirs})
endmacro()

# Create a target from find_package results.
#
# This is useful for packages with
# an 'old-way' find_package cmake routine.
# 
# Usage:
# 
# find_package(<name>)
# create_target(NAME <name> LIBRARIES <list of libs>  INCLUDE_DIR <list of includes>)
# target_link_libraries(some_other_target PRIVATE name)
# 
# Result : create a target <name>.
# some_other_target will be linked with <list of libs> and use <list of includes>
# to search for headers.
#
# See example in mechanics/CMakeLists.txt, for Bullet setup.
# 
function(create_target)
  set(oneValueArgs NAME)
  set(multiValueArgs LIBRARIES INCLUDE_DIRS)
  cmake_parse_arguments(target "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT TARGET ${target_NAME})
    add_library(${target_NAME} IMPORTED INTERFACE)
    set_property(TARGET ${target_NAME} PROPERTY INTERFACE_LINK_LIBRARIES ${target_LIBRARIES})
    if(target_INCLUDE_DIRS)
      set_target_properties(${target_NAME} PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${target_INCLUDE_DIRS}")
    endif()
  endif()
endfunction()


# Apply sanitizer options onto a given target
#
# Depends on user-defined variable USE_SANITIZER.
# 
# Might be:
# - asan : fast memory error detector
# - leaks : detects memory leaks
# - msan : detects uninitialized reads.
# - threads : detects data race
# - cfi : control flow integrity
# - undefined : detects the use of various features of C/C++ that are explicitly listed as resulting in undefined behaviour.
#
# Warning : do not combine options and notice that some of them may fail
# on MacOs ...
# 
# Ref :
# - http://www.stablecoder.ca/2018/10/30/full-cmake-helper-suite.html
# - https://clang.llvm.org/docs/AddressSanitizer.html
function(apply_sanitizer CURRENT_TARGET)

  unset(SANITIZER_OPTIONS)
  if(USE_SANITIZER MATCHES "asan") # Address sanitizer
    list(APPEND SANITIZER_OPTIONS "-fsanitize=address")
  elseif(USE_SANITIZER MATCHES "undefined") # Address sanitizer
    list(APPEND SANITIZER_OPTIONS "-fsanitize=undefined")
  elseif(USE_SANITIZER MATCHES "msan") # Memory sanitizer
    list(APPEND SANITIZER_OPTIONS "-fsanitize=memory")
    list(APPEND SANITIZER_OPTIONS "-fsanitize-memory-track-origins")
  elseif (USE_SANITIZER STREQUAL "thread") # Thread sanitizer (detect data race ...)
    list(APPEND SANITIZER_OPTIONS "-fsanitize=thread")
  elseif (USE_SANITIZER STREQUAL "Leak")
    list(APPEND SANITIZER_OPTIONS "-fsanitize=leak") # not implemented on macos (01/2020)
  elseif(USE_SANITIZER MATCHES "cfi") # control flow integrity
    list(APPEND SANITIZER_OPTIONS "-fsanitize=cfi")
    list(APPEND SANITIZER_OPTIONS "-flto")
    list(APPEND SANITIZER_OPTIONS "-B ${CLANG_LD_HACK}")
  endif()
  if(DEFINED SANITIZER_OPTIONS)
    list(APPEND SANITIZER_OPTIONS "-fno-omit-frame-pointer")
    message(STATUS "Activate sanitizer options (USE_SANITIZER=${USE_SANITIZER}) : ${SANITIZER_OPTIONS}")
  endif()

  
  target_compile_options(${CURRENT_TARGET} PUBLIC ${SANITIZER_OPTIONS})
  target_link_options(${CURRENT_TARGET} PUBLIC ${SANITIZER_OPTIONS})
  
endfunction()
