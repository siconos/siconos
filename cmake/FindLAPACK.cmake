
# - Find LAPACK library
# This module finds an installed library that implements the LAPACK 
# with its C interface (cblas) and the classical fortran interface.
# 
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use LAPACK
#  LAPACK_INCLUDE_DIRS - location of blas headers found by cmake
#  LAPACK_LIBRARY_DIR - location of the first blas lib found by cmake.
#  LAPACKE_HEADER - name of the header required for lapacke.
#  CLAPACK_HEADER - name of the header required for clapack. 
#  We first look for lapacke header (if it may be available in lapack implem) and then, if not found, for CLAPACK_HEADER.
#
# To find a specific blas set : 
#  WITH_LAPACK  if set checks only a specific implementation of blas which may be : mkl, openblas, atlas, accelerate, veclib, generic. 
# Or : 
#  LAPACK_DIR if set, cmake uses this path (absolute) to find a blas implementation in LAPACK_DIR/include and/or LAPACK_DIR/lib 
# 
# If neither WITH_LAPACK nor LAPACK_DIR is set, all implementation are searched in this order : mkl, openblas, accelerate, veclib, atlas, generic. 
# If available, pkg-config system is used to give hints to cmake search process.
#
if(NOT LAPACK_FOUND)
  
  include(CheckFunctionExists)
  include(CheckFortranFunctionExists)
  find_package(PkgConfig)


  # Check the language being used
  # If only C/C++, we check for cblas_... functions.
  # If Fortran and C, we check also the fortran interface.
  get_property( _LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES )
  if( _LANGUAGES_ MATCHES Fortran )
    set( _CHECK_FORTRAN TRUE )
  elseif( (_LANGUAGES_ MATCHES C) OR (_LANGUAGES_ MATCHES CXX) )
    set( _CHECK_FORTRAN FALSE )
  else()
    if(LAPACK_FIND_REQUIRED)
      message(FATAL_ERROR "FindLAPACK requires Fortran, C, or C++ to be enabled.")
    else()
      message(STATUS "Looking for LAPACK... - NOT found (Unsupported languages)")
      return()
    endif()
  endif()
  
  macro(check_lapack_Libraries LIBRARIES _prefix _name _flags _list _blas _threads)
    # This macro checks for the existence of the combination of fortran libraries
    # given by _list.  If the combination is found, this macro checks (using the
    # Check_Fortran_Function_Exists macro) whether can link against that library
    # combination using the name of a routine given by _name using the linker
    # flags given by _flags.  If the combination of libraries is found and passes
    # the link test, LIBRARIES is set to the list of complete library paths that
    # have been found.  Otherwise, LIBRARIES is set to FALSE.
    
    # N.B. _prefix is the prefix applied to the names of all cached variables that
    # are generated internally and marked advanced by this macro.
    set(_libdir ${ARGN})
    
    set(_libraries_work TRUE)
    set(${LIBRARIES} "")
    set(_combined_name)
    
    ## If no extra argument was given to the macro, default search path is
    ## filled with environment variables.
    if (NOT _libdir)
      if (WIN32)
	set(_libdir LIB)
      elseif (APPLE)
	set(_libdir ENV DYLD_LIBRARY_PATH)
      else ()
	set(_libdir ENV LD_LIBRARY_PATH)
      endif ()
    endif ()
    foreach(_library ${_list})
      set(_combined_name ${_combined_name}_${_library})
      
      if(_libraries_work)
	
	# HINTS are checked before PATHS, that's why we call
	# find_library twice, to give priority to LD_LIBRARY_PATH or user-defined paths
	# over pkg-config process.
	# This first call should find _library in env. variables. 
	find_library(${_prefix}_${_library}_LIBRARIES
	  NAMES ${_library}
	  PATHS ${_libdir}
          PATH_SUFFIXES atlas
	  NO_DEFAULT_PATH
	  )
	find_library(${_prefix}_${_library}_LIBRARIES
	  NAMES ${_library}
          PATH_SUFFIXES atlas
	  )
	## If search fails, we try with pkg-config
	if(NOT ${_prefix}_${_library}_LIBRARIES)
	  set(${_prefix}_${_library}_LIBRARIES "")
	  message(STATUS "Try to find ${_library} using pkg-config")
	  pkg_check_modules(PC_${_library} QUIET ${_library})
	  foreach(PC_LIB ${PC_${_library}_LIBRARIES})
	    find_library(${PC_LIB}_LIBRARY
	      NAMES ${PC_LIB}
	      HINTS ${PC_${_library}_LIBRARY_DIRS} 
	      )
	    if (NOT ${PC_LIB}_LIBRARY)
	      message(FATAL_ERROR "Something is wrong in your pkg-config file - lib ${PC_LIB} not found in ${PC_${_library}_LIBRARY_DIRS}")
	    endif (NOT ${PC_LIB}_LIBRARY)
	    list(APPEND ${_prefix}_${_library}_LIBRARIES ${${PC_LIB}_LIBRARY}) 
	  endforeach(PC_LIB)
	  ## pkg-config may give some hints about headers location
	  set(INCLUDE_DIR_HINTS ${PC_${_library}_INCLUDE_DIRS})
	endif()
	mark_as_advanced(${_prefix}_${_library}_LIBRARIES)
	set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARIES})
	set(_libraries_work ${${_prefix}_${_library}_LIBRARIES})
      endif()
    endforeach()
    if(_libraries_work)
      # Test this combination of libraries.
      if(UNIX AND BLA_STATIC)
	set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
      else()
	set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
      endif()
      # add gfortran if we have static libs + gfortran
      #      if (BLA_STATIC AND CMAKE_COMPILER_IS_GNUG77)
      #   if (NOT GFORTRAN_LIB)
      #     set(GFORTRAN_LIB "gfortran")
      #   endif()
      #   set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GFORTRAN_LIB})
      #  endif()
      ## First we check c interface
      if(${_prefix} STREQUAL "LAPACKE")
        check_function_exists("LAPACKE_${_name}" ${_prefix}${_combined_name}_WORKS)
      else()
	check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
      endif()
      # and then, if required, fortran interface
      if (_CHECK_FORTRAN)
	check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS_F)
	if(NOT ${${_prefix}${_combined_name}_WORKS_F})
	  set(${_prefix}${_combined_name}_WORKS FALSE)
	endif()
      endif()
      set(CMAKE_REQUIRED_LIBRARIES)
      mark_as_advanced(${_prefix}${_combined_name}_WORKS)
      set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
    endif()
    if(_libraries_work)
      set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
    else()
      set(${LIBRARIES} FALSE)
    endif()
  endmacro()
  
  if (BLA_STATIC)
    if (WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif ()
    if (APPLE)
      set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else ()
      set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif ()
  else ()
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # for ubuntu's libblas3gf and liblapack3gf packages
      set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
    endif ()
  endif ()

  ## First of all, we need blas ## 
  if(NOT BLAS_FOUND)
    if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
      message(STATUS "Warning : BLAS not found on your system. Lapack will not be searched.")
    else()
      message(FATAL_ERROR "You need to find BLAS before any attempt to find LAPACK.")
    endif()
  endif()
  
  if(BLAS_FOUND)
    #### Start Lapack search process ####
    # I'm not sure what this is supposed to achieve, but it prevents us from using an hint on the command line ... -- xhub
    #    set(WITH_LAPACK "" CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/generic]")
    set(LAPACK_DIR "" CACHE PATH "lapack implementation location.")
    
    ## We first check the blas implementation : if it is mkl or accelerate
    ## then lapack must be in mkl/accelerate. 
    ## 
    # MKL Intel  - Everything is already determined during blas search. 
    if(WITH_BLAS STREQUAL "mkl")
      set(WITH_LAPACK "mkl" CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE)
      set(LAPACKE_HEADER mkl_lapacke.h)
      set(HAS_LAPACKE TRUE)
      set(LAPACK_VERSION ${MKL_VERSION})
      set(LAPACK_LIBRARY_DIR ${BLAS_LIBRARY_DIR} CACHE PATH "Lapack libraries location." FORCE)
      set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS} CACHE PATH "Lapack header location." FORCE)
      set(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES})
      ## Apple accelerate 
    elseif(WITH_BLAS STREQUAL "accelerate")
      set(WITH_LAPACK "accelerate" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE)
      set(CLAPACK_HEADER clapack.h)
      set(LAPACK_LIBRARY_DIR ${BLAS_LIBRARY_DIR} CACHE PATH "Lapack libraries location." FORCE)
      
    else() ## Other blas implementations
      # if nothing set by user for lapack implementation, we copy blas one's. 
      if((NOT WITH_LAPACK) AND (NOT LAPACK_DIR))
	set(WITH_LAPACK ${WITH_BLAS} CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE )

	# if LAPACK_DIR is provided, we add it to CMAKE_PREFIX_PATH. This should be the first place to be 
	# searched for lapack libraries and headers.
      elseif((NOT WITH_LAPACK) AND LAPACK_DIR)
	set(CMAKE_PREFIX_PATH ${LAPACK_DIR})
	string(TOLOWER ${LAPACK_DIR} lowerlapdir)
	## We use LAPACK_DIR to set WITH_LAPACK, if possible. 
	if(lowerlapdir MATCHES "atlas.*") 
	  set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	elseif(lowerlapdir MATCHES "openblas.*")
	  set(WITH_LAPACK "openblas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	endif()
      endif()
      
    endif()
    
    ## Now we search and check lapack libraries and headers. 
    ## MKL : already done thanks to blas
    
    ## OpenBLAS ##
    # this can come in 3 flavors: no LAPACK, fortran LAPACK, LAPACKE
    # the problem with the second option is that we have no headers ...
    # the last case occurs for instance on Fedora (up to 21 at time of writing)
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "openblas")))
      message(STATUS "Try to find lapack in openblas ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACKE
	cheev
	""
	"openblas"
	"${BLAS_LIBRARIES}"
	"")
      if(LAPACK_LIBRARIES)
	set(WITH_LAPACK "openblas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(LAPACKE_HEADER lapacke.h)
	set(LAPACK_INCLUDE_SUFFIXES "openblas;openblas/lapacke")
      else()
	set(WITH_LAPACK "" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      endif(LAPACK_LIBRARIES)
    endif()

    ## Apple Framework ## 
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "accelerate")))
      message(STATUS "Try to find lapack in Accelerate framework ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"Accelerate"
	"${BLAS_LIBRARIES}"
  	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "accelerate" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
	set(LAPACK_INCLUDE_SUFFIXES Headers Frameworks)
      endif (LAPACK_LIBRARIES)
    endif()
    
    ## Generic LAPACKE library ##
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "lapacke") OR (WITH_LAPACK STREQUAL "generic")))
      message(STATUS "Try to find a generic lapacke ...")
      if(LAPACKE_NAME)
	message(STATUS "Using user-provided name ${LAPACKE_NAME}")
      else(LAPACKE_NAME)
	set(LAPACKE_NAME "lapacke")
      endif(LAPACKE_NAME)
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACKE
	cheev
	""
	"${LAPACKE_NAME};lapack"
	"${BLAS_LIBRARIES}"
	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "generic" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(LAPACKE_HEADER lapacke.h)
      endif()
    endif()

    ## Atlas ## 
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "atlas")))
      message(STATUS "Try to find lapack in atlas ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"lapack_atlas;lapack"
	"${BLAS_LIBRARIES}"
	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
	set(LAPACKE_HEADER FALSE)
	set(LAPACK_INCLUDE_SUFFIXES atlas)
      else()
	set(WITH_LAPACK "" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      endif()
    endif()

    ## Generic LAPACK library ##
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "generic")))
      message(STATUS "Try to find a generic lapack ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"lapack"
	"${BLAS_LIBRARIES}"
	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "generic" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
	set(LAPACK_INCLUDE_SUFFIXES atlas) # For debian or ubuntu ...
      endif (LAPACK_LIBRARIES)
    endif()
    
    if(LAPACK_LIBRARIES)
      set(LAPACK_FOUND TRUE)
    else()
      set(LAPACK_FOUND FALSE)
    endif()
    
    ## Now the headers ...
    if(LAPACK_FOUND)
      unset(LAPACK_FOUND CACHE)
      message(STATUS "Lapack libraries have been found : ${LAPACK_LIBRARIES}. We now turn to headers.")

      if(NOT LAPACK_LIBRARY_DIR)
	# Get an hint for the location of lapack headers : probably LAPACK_LIBRARIES/../include or LAPACK_DIR if given
	list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
	get_filename_component(_bdir ${LAPACK_LIB} PATH)
	set(LAPACK_LIBRARY_DIR ${_bdir} CACHE PATH "Lapack libraries location." FORCE)
      endif()
      
      # Either LAPACK_DIR has been set by user or we use LAPACK_LIBRARIES location to set it. 
      if(NOT LAPACK_DIR)
	get_filename_component(_bdir ${LAPACK_LIBRARY_DIR} PATH)
	set(LAPACK_DIR ${_bdir} CACHE PATH "Lapack implementation location." FORCE)
      endif()
      
      set(LAPACK_INC_DIR ${LAPACK_DIR}/include)
      ## Use LAPACK_DIR to set path for cmake search. 
      ## If pkg-config has been used, maybe it gives some hints about headers location ...
      set(CMAKE_INCLUDE_PATH ${LAPACK_INC_DIR} ${INCLUDE_DIR_HINTS})
      set(CMAKE_PREFIX_PATH ${LAPACK_DIR})

      ## Note Franck : it seems that find_path process does not work as expected on Macosx : it does not respect the search sequence described
      # in cmake doc (i.e. CMAKE_PREFIX, then HINTS, PATHS ... ) and always turn to find apple framework cblas if
      # NO_DEFAULT_PATH is not set. 
      if(APPLE) # First check in HINTS, no default, then global search.
	if(LAPACKE_HEADER)
	  find_path(LAPACK_INCLUDE_DIRS
	    NAMES ${LAPACKE_HEADER}
	    HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    NO_DEFAULT_PATH
	    )
	  find_path(LAPACK_INCLUDE_DIRS
	    NAMES ${LAPACKE_HEADER}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_LAPACKE 1 CACHE BOOL "Lapacke interface is available.")
	  endif()

	endif() # if lapacke.h not needed or not found search for clapack.h ...
	if(NOT LAPACK_INCLUDE_DIRS)
	  find_path(LAPACK_INCLUDE_DIRS
	    NAMES ${CLAPACK_HEADER}
	    HINTS ${LAPACK_LIBRARIES}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    NO_DEFAULT_PATH
	    )
	  find_path(LAPACK_INCLUDE_DIRS
	    NAMES ${CLAPACK_HEADER}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
	  else()
	    if(LAPACK_FIND_REQUIRED)
	      message(FATAL_ERROR "Lapack headers not found.")
	    endif()
	  endif()
	endif()

      else(APPLE) # The case which is supposed to always work
	if(LAPACKE_HEADER)
	  find_path(LAPACK_INCLUDE_DIRS 
	    NAMES ${LAPACKE_HEADER}
	    HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_LAPACKE 1 CACHE BOOL "Lapacke interface is available.")
	  endif()
	endif()
	# if lapacke.h not needed or not found search for clapack.h ...
	if(NOT LAPACK_INCLUDE_DIRS)
	  find_path(LAPACK_INCLUDE_DIRS 
	    NAMES ${CLAPACK_HEADER}
	    HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
	  endif()

	endif()
      endif()
      if(LAPACK_INCLUDE_DIRS)
	set(LAPACK_FOUND 1 CACHE BOOL "Found lapack headers dir.")
      endif()
      
    endif()
  endif()

   if(LAPACK_FOUND) # SiconosConfig.h setup
     if(NOT LAPACK_LIBRARIES MATCHES "lapacke.*") # atlas can be the BLAS provider
      ## If lapack was found as "generic" but is part of atlas. 
      if(LAPACK_LIBRARIES MATCHES "atlas.*") 
         set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
       endif()
     endif()

    if(WITH_LAPACK STREQUAL "mkl")
      set(HAS_MKL_LAPACKE 1 CACHE BOOL "Blas comes from Intel MKL.")
      set(LAPACK_SUFFIX)
      set(LAPACK_PREFIX "LAPACKE_")

    elseif(WITH_LAPACK STREQUAL "accelerate")
      set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
      set(LAPACK_SUFFIX "_")
      set(LAPACK_PREFIX)

    elseif(WITH_LAPACK STREQUAL "atlas")
      set(HAS_ATLAS_LAPACK 1 CACHE BOOL "Blas  comes from Atlas framework ")
      set(LAPACK_SUFFIX)
      set(LAPACK_PREFIX "clapack_")
      set(HAS_LAPACKE FALSE)

    elseif(WITH_LAPACK STREQUAL "openblas")
      set(HAS_OpenBLAS_LAPACK 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")
      set(LAPACK_SUFFIX)
      if(HAS_LAPACKE)
        set(LAPACK_PREFIX "LAPACKE_")
      else()
      endif()
    else()
      set(HAS_GenericCLAPACK 1 CACHE BOOL "Lapack is available from an unknown version.")
      if(HAS_LAPACKE)
	set(LAPACK_PREFIX "LAPACKE_")
      else()
	set(LAPACK_PREFIX "clapack_")
      endif()
    endif()
  endif()

  if(HAS_LAPACKE)
    set(LAPACK_HEADER ${LAPACKE_HEADER} CACHE STRING "Lapack header name.")
  else()
    set(LAPACK_HEADER ${CLAPACK_HEADER} CACHE STRING "Lapack header name.")
  endif()    
    # === Extra Checks for lapack functions ===
  # We check only the functions that are known to be un-implemented
  # in some lapack libraries (namely atlas ...)
  # This is probably a temporary check since it's likely
  # we will stop atlas checking for lapack? 
  ## Check if FuncName is available in lapack lib (i.e in one of LAPACK_LIBRARIES)
  ## and if FuncName symbol is present in file Header. 
  # If both are true, result is true.
  function(check_clapack_has_function genericName FuncName Header result)

    check_function_exists(${FuncName} TEST_FUNC)
    check_symbol_exists(${FuncName} ${Header} TEST_HEAD)

    if(TEST_HEAD AND TEST_FUNC)
      set(${result} 1 CACHE BOOL "${genericName} is available in lapack.")
    endif()

    unset(TEST_HEAD)
    unset(TEST_FUNC)
  endfunction()


  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
  if (BLA_STATIC AND CMAKE_COMPILER_IS_GNUG77 AND NOT MSVC)
    if (NOT GFORTRAN_LIB)
      set(GFORTRAN_LIB "gfortran")
    endif()
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GFORTRAN_LIB})
  endif()

  if(MSVC AND HAS_LAPACKE)
    SET(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS} "-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE")
    SET(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${EXTRA_LAPACK_LIB})
  endif(MSVC AND HAS_LAPACKE)

  ## dgesvd ##
  unset(HAS_LAPACK_DGESVD)
  set(GENERIC_NAME "DGESVD")
  set(FUNC_NAME "${LAPACK_PREFIX}dgesvd${LAPACK_SUFFIX}")
  check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGESVD)

  ## dgels ##
  unset(HAS_LAPACK_DGELS)
  set(GENERIC_NAME "DGELS")
  set(FUNC_NAME "${LAPACK_PREFIX}dgels${LAPACK_SUFFIX}")
  check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGELS)

  ## dtrtrs ##
  unset(HAS_LAPACK_DTRTRS)
  set(GENERIC_NAME "DTRTRS")
  set(FUNC_NAME "${LAPACK_PREFIX}dtrtrs${LAPACK_SUFFIX}")
  check_clapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DTRTRS)

  ## Final settings ...
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "Lapack libraries.")

  if(NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "Cannot find a library with LAPACK API. Please specify library location using LAPACK_DIR option or set your environment variables properly.")
  endif (NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
  if(NOT LAPACK_FIND_QUIETLY)
    if(LAPACK_FOUND)
      message(STATUS "Found a library with LAPACK API (${WITH_LAPACK}).")
      message("    Lapack libraries : ${LAPACK_LIBRARIES}.")
      message("    Lapack headers : ${LAPACK_HEADER} in ${LAPACK_INCLUDE_DIRS}.")
    else(LAPACK_FOUND)
      message(STATUS "Cannot find a library with LAPACK API. Maybe you can try again using LAPACK_DIR option or set your environment variables properly.")
    endif(LAPACK_FOUND)
  endif(NOT LAPACK_FIND_QUIETLY)
endif()




