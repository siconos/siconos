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
    if(BLAS_FIND_REQUIRED)
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
    set(${LIBRARIES})
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
	
	# HINTS are checked before PATHS, that's why we call
	# find_library twice, to give priority to LD_LIBRARY_PATH or user-defined paths
	# over pkg-config process.
	find_library(${_prefix}_${_library}_LIBRARY
	  NAMES ${_library}
	  PATHS ${_libdir}
	  NO_DEFAULT_PATH
	  )
	pkg_check_modules(PC_LIBRARY QUIET ${_library})
	PRINT_VAR(PC_LIBRARY_LIBDIR)
	PRINT_VAR(PC_LIBRARY_LIBRARY_DIRS)
	find_library(${_prefix}_${_library}_LIBRARY
	  NAMES ${_library}
	  HINTS ${PC_LIBRARY_LIBDIR} ${PC_LIBRARY_LIBRARY_DIRS} 
	  )
	PRINT_VAR(${_prefix}_${_library}_LIBRARY)
	mark_as_advanced(${_prefix}_${_library}_LIBRARY)
	set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
	set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
      endif()
    endforeach()
    if(_libraries_work)
      # Test this combination of libraries.
      if(UNIX AND BLA_STATIC)
	set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
      else()
	set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
      endif()
      ## First we check c interface
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
      # and then, if required, fortran interface
      if (_CHECK_FORTRAN)
	check_fortran_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS_F)
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
  
  ## First of all, we need blas ## 
  if(NOT BLAS_FOUND)
    if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
      compile_with(BLAS)
    else()
      compile_with(BLAS REQUIRED)
    endif()
  endif()
  
  if(BLAS_FOUND)
    #### Start Blas search process ####
    set(WITH_LAPACK "" CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]")
    set(LAPACK_DIR "" CACHE PATH "lapack implementation location.")
    
    ## We first check the blas implementation ##
    # MKL Intel  - Everything is already determined during blas search. 
    if(WITH_BLAS STREQUAL "mkl")
      set(WITH_LAPACK "mkl" CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
      set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE)
      set(LAPACKE_HEADER mkl_lapacke.h)
      set(LAPACK_LIBRARIES ${MKL_LAPACK_LIBRARIES})
      set(LAPACK_INCLUDE_DIR ${MKL_INCLUDE_DIR})
      set(LAPACK_VERSION ${MKL_VERSION})
      set(LAPACK_LIBRARY_DIR ${BLAS_LIBRARY_DIR} CACHE PATH "Lapack libraries location." FORCE)
      
     ## Appel accelerate or veclib 
    elseif(WITH_BLAS STREQUAL "accelerate")
      set(WITH_LAPACK "accelerate" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
      set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE)
      set(CLAPACK_HEADER clapack.h)
      set(LAPACK_LIBRARY_DIR ${BLAS_LIBRARY_DIR} CACHE PATH "Lapack libraries location." FORCE)

    elseif(WITH_BLAS STREQUAL "veclib")
      set(WITH_LAPACK "veclib" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
      set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE)
      set(CLAPACK_HEADER clapack.h)
      set(LAPACK_LIBRARY_DIR ${BLAS_LIBRARY_DIR} CACHE PATH "Lapack libraries location." FORCE)
      
    else() ## Other blas implementations
      # if nothing set by user for lapack implementation, we copy blas one's. 
      if((NOT WITH_LAPACK) AND (NOT LAPACK_DIR))
	set(WITH_LAPACK ${WITH_BLAS} CACHE STRING "Lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(LAPACK_DIR ${BLAS_DIR} CACHE PATH "lapack implementation location." FORCE )
      elseif((NOT WITH_LAPACK) AND LAPACK_DIR)
	set(CMAKE_PREFIX_PATH ${LAPACK_DIR})
	string(TOLOWER ${LAPACK_DIR} lowerlapdir)
	if(lowerblasdir MATCHES "atlas.*") 
	  set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	elseif(lowerblasdir MATCHES "openblas.*")
	  set(WITH_LAPACK "openblas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	endif()
      endif()
	
    endif()

    ## Now we search and check lapack libraries and headers. 
    ## MKL : already done thanks to blas
    
    ## OpenBLAS ##
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "openblas")))
      message(STATUS "Try to find lapack in openblas ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"openblas"
	"${BLAS_LIBRARIES}"
	"")
      if(LAPACK_LIBRARIES)
	set(WITH_LAPACK "openblas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(LAPACKE_HEADER lapacke.h)
	set(CLAPACK_HEADER clapack.h)
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
	set(WITH_LAPACK "accelerate" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
      endif (LAPACK_LIBRARIES)
    endif()
    if((NOT LAPACK_LIBRARIES)
	AND ((NOT WITH_LAPACK) OR (WITH_LAPACK STREQUAL "veclib")))
      message(STATUS "Try to find lapack in VecLib framework ...")
      check_lapack_libraries(
	LAPACK_LIBRARIES
	LAPACK
	cheev
	""
	"vecLib"
	"${BLAS_LIBRARIES}"
	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "veclib" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
      endif (LAPACK_LIBRARIES)
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
	"lapack"
	"${BLAS_LIBRARIES}"
	"")
      if (LAPACK_LIBRARIES)
	set(WITH_LAPACK "atlas" CACHE STRING "lapack implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(CLAPACK_HEADER clapack.h)
	set(LAPACK_INCLUDE_SUFFIXES atlas)
      endif (LAPACK_LIBRARIES)
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
	set(WITH_LAPACK "generic" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/veclib/generic]" FORCE)
	set(LAPACKE_HEADER lapacke.h)
	set(CLAPACK_HEADER clapack.h)
	set(LAPACK_INCLUDE_SUFFIXES atlas) # For debian or ubuntu ...
      endif (LAPACK_LIBRARIES)
    endif()
    
    if(LAPACK_LIBRARIES)
      set(LAPACK_FOUND TRUE)
    else(LAPACK_LIBRARIES)
      set(LAPACK_FOUND FALSE)
    endif(LAPACK_LIBRARIES)
    
    ## Now the headers ...
    if(LAPACK_FOUND)
      set(LAPACK_FOUND 0)
      if(NOT LAPACK_LIBRARY_DIR)
	# Get an hint for the location of lapack headers : probably LAPACK_LIBRARIES/../include or LAPACK_DIR if given
	list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
	get_filename_component(_bdir ${LAPACK_LIB} PATH)
	set(LAPACK_LIBRARY_DIR ${_bdir} CACHE PATH "Lapack libraries location." FORCE)
      endif()
      
      if(NOT LAPACK_DIR)
	get_filename_component(_bdir ${LAPACK_LIBRARY_DIR} PATH)
	set(LAPACK_DIR ${_bdir} CACHE PATH "Lapack implementation location." FORCE)
      endif()
      
      PRINT_VAR(LAPACK_LIBRARY_DIR)
      PRINT_VAR(LAPACK_DIR)
       
      set(LAPACK_INC_DIR ${LAPACK_DIR}/include)
      set(CMAKE_INCLUDE_PATH ${LAPACK_INC_DIR})
      set(CMAKE_PREFIX_PATH ${LAPACK_DIR})
      PRINT_VAR(LAPACK_INC_DIR)

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
	    HINTS ${LAPACK_DIR} ${LAPACK_INC_DIR}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    NO_DEFAULT_PATH
	    )
	  find_path(LAPACK_INCLUDE_DIRS
	    NAMES ${CLAPACK_HEADER}
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
	  endif()
	endif()

	#set(LAPACK_INCLUDE_DIRS ${_dir} CACHE STRING "Lapack header location." FORCE)
	
      else() # The case which is supposed to always work
	if(LAPACKE_HEADER)
	  find_path(LAPACK_INCLUDE_DIRS 
	    NAMES ${LAPACKE_HEADER}
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
	    PATH_SUFFIXES ${LAPACK_INCLUDE_SUFFIXES}
	    )
	  if(LAPACK_INCLUDE_DIRS)
	    set(HAS_CLAPACK 1 CACHE BOOL "clapack interface is available.")
	  endif()

	endif()
      endif()
      PRINT_VAR(LAPACK_DIR)
      PRINT_VAR(LAPACK_INCLUDE_DIRS)
      if(LAPACK_INCLUDE_DIRS)
	set(LAPACK_FOUND 1)
      endif()
      
    endif()
  endif()

   if(LAPACK_FOUND) # NumericsConfig.h setup
     
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

    elseif(WITH_LAPACK STREQUAL "openblas")
      set(HAS_OpenBLAS_LAPACK 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")
      set(LAPACK_SUFFIX)
      if(HAS_LAPACKE)
	set(LAPACK_PREFIX "LAPACKE_")
      else()
	set(LAPACK_PREFIX "clapack_")
      endif()
    else()
      set(HAS_GenericCLAPACK 1 CACHE BOOL "Blas is available from an unknown version.")
      if(HAS_LAPACKE)
	set(LAPACK_PREFIX "LAPACKE_")
      else()
	set(LAPACK_PREFIX "clapack_")
      endif()
    endif()
  endif()

  if(NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "Cannot find a library with LAPACK API. Please specify library location using LAPACK_DIR option or set your environment variables properly.")
  endif (NOT LAPACK_FOUND AND LAPACK_FIND_REQUIRED)
  if(NOT LAPACK_FIND_QUIETLY)
    if(LAPACK_FOUND)
      message(STATUS "Found a library with LAPACK API (${BLA_VENDOR}).")
    else(LAPACK_FOUND)
      message(STATUS "Cannot find a library with LAPACK API. Maybe you can try again using LAPACK_DIR option or set your environment variables properly.")
    endif(LAPACK_FOUND)
  endif(NOT LAPACK_FIND_QUIETLY)
  
  # === Extra Checks for lapack functions ===
  # We check only the functions that are known to be un-implemented
  # in some lapack libraries (namely atlas ...)
  # This is probably a temporary check since it's likely
  # we will stop atlas checking for lapack? 
  ## Check if FuncName is available in lapack lib (i.e in one of LAPACK_LIBRARIES)
  ## and if FuncName symbol is present in file Header. 
  # If both are true, result is true.
  function(check_lapack_has_function genericName FuncName Header result)
    
    check_function_exists(${FuncName} TEST_FUNC)
    check_symbol_exists(${FuncName} ${Header} TEST_HEAD)
    
    if(TEST_HEAD AND TEST_FUNC)
      set(${result} 1 CACHE BOOL "${genericName} is available in lapack.")
    endif()
    
  endfunction()
  
  
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
  ## dgesvd ##
  unset(HAS_LAPACK_DGESVD)
  set(GENERIC_NAME "DGESVD")
  set(FUNC_NAME "${LAPACK_PREFIX}dgesvd${LAPACK_SUFFIX}")
  if(HAS_LAPACKE)
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACKE_HEADER} HAS_LAPACK_DGESVD)
  else()
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${CLAPACK_HEADER} HAS_LAPACK_DGESVD)
  endif()
  ## dgels ##
  unset(HAS_LAPACK_DGELS)
  set(GENERIC_NAME "DGELS")
  set(FUNC_NAME "${LAPACK_PREFIX}dgels${LAPACK_SUFFIX}")
  if(HAS_LAPACKE)
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACKE_HEADER} HAS_LAPACK_DGELS)
  else()
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${CLAPACK_HEADER} HAS_LAPACK_DGELS)
  endif()

  ## dtrtrs ##
  unset(HAS_LAPACK_DTRTRS)
  set(GENERIC_NAME "DTRTRS")
  set(FUNC_NAME "${LAPACK_PREFIX}dtrtrs${LAPACK_SUFFIX}")
  if(HAS_LAPACKE)
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACKE_HEADER} HAS_LAPACK_DTRTRS)
  else()
    check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${CLAPACK_HEADER} HAS_LAPACK_DTRTRS)
  endif()





endif()




