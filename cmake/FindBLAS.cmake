# - Find BLAS library
# This module finds an installed library that implements the BLAS
# with its C interface (cblas) and the classical fortran interface.
#
#
# This module sets the following variables:
#  BLAS_FOUND - set to true if a library implementing the BLAS interface
#    is found
#  BLAS_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use BLAS
#  BLAS_INCLUDE_DIRS - location of blas headers found by cmake
#  BLAS_LIBRARY_DIR - location of the first blas lib found by cmake.
#  BLAS_HEADER - name(s) of the header(s) required for cblas.
#
# To find a specific blas set :
#  WITH_BLAS  if set checks only a specific implementation of blas which may be : mkl, openblas, atlas, accelerate, veclib, generic.
# Or :
#  BLAS_DIR if set, cmake uses this path (absolute) to find a blas implementation in BLAS_DIR/include and/or BLAS_DIR/lib
#
# If neither WITH_BLAS nor BLAS_DIR is set, all implementation are searched in this order : mkl, openblas, accelerate, veclib, atlas, generic.
# If available, pkg-config system is used to give hints to cmake search process.
#


if(NOT BLAS_FOUND)

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
      message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
    else()
      message(STATUS "Looking for BLAS... - NOT found (Unsupported languages)")
      return()
    endif()
  endif()

  macro(check_blas_libraries LIBRARIES _prefix _name _flags _list _thread)
    # This macro checks for the existence of the combination of libraries
    # given by _list.  If the combination is found, this macro checks (using the
    # Check_Fortran_Function_Exists and/or Check_function_exists macros) whether
    # we can link against that library
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
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_thread})
      # add gfortran if we have static libs + gfortran
      #      if (BLA_STATIC AND CMAKE_COMPILER_IS_GNUG77)
      #        if (NOT GFORTRAN_LIB)
      #   set(GFORTRAN_LIB "gfortran")
      #  endif()
      #  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${GFORTRAN_LIB})
      #endif()
      #else()
      ## First we check cblas interface
      if (_libraries_work MATCHES "mwblas")
        set(_symbol_name "${_name}_")
      else()
        set(_symbol_name "cblas_${_name}")
      endif()
      check_function_exists("${_symbol_name}" ${_prefix}${_combined_name}_WORKS)
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
    if(NOT _libraries_work)
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

  #### Start Blas search process ####
  set(WITH_BLAS "" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]")
  set(BLAS_DIR "" CACHE PATH "Blas implementation location.")

  ## If BLAS_DIR is provided by user, we add it to CMAKE_PREFIX_PATH. This should be the first place to be
  # searched for blas libraries and headers.
  if(BLAS_DIR)
    set(CMAKE_PREFIX_PATH ${BLAS_DIR})
    string(TOLOWER ${BLAS_DIR} lowerblasdir)
  else()
    set(lowerblasdir)
  endif()

  ## Then we use BLAS_DIR to get hints on the blas implementation and set WITH_BLAS
  if(BLAS_DIR MATCHES "mkl.*")
    set(WITH_BLAS "mkl" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
    set(INTEL_MKL_DIR ${BLAS_DIR})
    set(INTEL_COMPILER_DIR "${BLAS_DIR}/..")
  elseif(lowerblasdir MATCHES "atlas.*")
    set(WITH_BLAS "atlas" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
  elseif(lowerblasdir MATCHES "openblas.*")
    set(WITH_BLAS "openblas" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
  elseif(lowerblasdir MATCHES "accelerate.framework.*")
    set(WITH_BLAS "accelerate" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
  elseif(lowerblasdir MATCHES "veclib.framework.*")
    set(WITH_BLAS "veclib" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
  endif()

  message(STATUS "Start blas search with BLAS_DIR =${BLAS_DIR} and WITH_BLAS=${WITH_BLAS}.")


  ## Now we will search for blas libraries in different implementations.
  ## If WITH_BLAS is not set, the first found will be set.
  ## We check for the fortran AND the c (cblas_) interface.
  ## If something is found, we MUST set:
  # BLAS_LIBRARIES : list of libraries to link with
  # WITH_BLAS : name of the implementation found (mkl/openblas/atlas/accelerate/generic)
  # BLAS_HEADER : name of the header for cblas interface.

  set(BLAS_LIBRARIES)
  ## Intel MKL ##
  if((NOT BLAS_LIBRARIES)
      AND ((NOT WITH_BLAS) OR (WITH_BLAS STREQUAL "mkl")))
    message(STATUS "Try to find blas in intel/mkl ...")
    find_package(MKL)
    if(MKL_FOUND)
      set(WITH_BLAS "mkl" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_LIBRARIES ${MKL_LIBRARIES})
      ## FindMKL also sets MKL_LAPACK_LIBRARIES
      set(INCLUDE_DIR_HINTS ${MKL_INCLUDE_DIR})
      set(BLAS_VERSION ${MKL_VERSION})
      set(BLAS_HEADER mkl_cblas.h CACHE STRING "Blas header name")
    endif(MKL_FOUND)
  endif()

  ## OpenBLAS ##
  if((NOT BLAS_LIBRARIES)
    AND ((NOT WITH_BLAS) OR (WITH_BLAS MATCHES "openblas")))
    if(WITH_BLAS MATCHES "openblas")
      set(OPENBLAS_LIB_NAME ${WITH_BLAS})
    else()
      set(OPENBLAS_LIB_NAME "openblas")
    endif()
    message(STATUS "Try to find blas in openblas ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "${OPENBLAS_LIB_NAME}"
      "")
    if(BLAS_LIBRARIES)
      set(WITH_BLAS "openblas" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER cblas.h CACHE STRING "Blas header name")
      # calling FORTRAN complex function from C++ is not easy, see
      # https://software.intel.com/sites/products/documentation/hpc/mkl/mkl_userguide_lnx/GUID-A0908E50-19D7-44C1-A068-44036B466BC7.htm
      # http://www.netlib.org/blas/blast-forum/cinterface.pdf
      # to ensure that std::complex<double> can be casted to double* and produces good results we should enfore std=c++11
      # Right now openblas ships with a cblas.h header which uses double* as type for pointer to double complex, which makes ublas cry since it
      # want to pass std::complex<double>*
      #set(BLAS_INCLUDE_SUFFIXES cblas_header)
      #set(INCLUDE_DIR_HINTS ${CMAKE_SOURCE_DIR}/externals/blas_lapack)

      get_filename_component(_bdir ${BLAS_LIBRARIES} DIRECTORY)
      set(INCLUDE_DIR_HINTS "${_bdir}/../include/openblas/")
      set(OPENBLAS_FOUND 1)

    endif(BLAS_LIBRARIES)
  endif()

  ## Apple Framework ##
  if((NOT BLAS_LIBRARIES)
      AND ((NOT WITH_BLAS) OR (WITH_BLAS STREQUAL "accelerate")))
    message(STATUS "Try to find blas in Accelerate framework ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "Accelerate"
      "")
    if (BLAS_LIBRARIES)
      set(WITH_BLAS "accelerate" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER Accelerate.h cblas.h CACHE STRING "Blas header name(s)")
      set(BLAS_INCLUDE_SUFFIXES Headers Frameworks)
    endif (BLAS_LIBRARIES)
  endif()

  ## Atlas ##
  ## Rq : atlas is the "default" cblas environmment in some linux distrib (Debian, Ubuntu, Fedora ...)
  if((NOT BLAS_LIBRARIES)
      AND ((NOT WITH_BLAS) OR (WITH_BLAS STREQUAL "atlas")))
    message(STATUS "Try to find blas in atlas ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "cblas;atlas;f77blas"
      "")
    if (BLAS_LIBRARIES)
      set(WITH_BLAS "atlas" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER cblas.h CACHE STRING "Blas header name")
      set(BLAS_INCLUDE_SUFFIXES atlas)
    endif (BLAS_LIBRARIES)
  endif()

  ## MATLAB ##
  if((NOT BLAS_LIBRARIES)
    AND ((NOT WITH_BLAS) OR (WITH_BLAS MATCHES "MATLAB")))
    set(MATLAB_LIB_NAME "mwblas")
    message(STATUS "Try to find blas in MATLAB ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      dgemm
      ""
      "${MATLAB_LIB_NAME}"
      "")
    if(BLAS_LIBRARIES)
      set(WITH_BLAS "MATLAB" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER blas.h CACHE STRING "Blas header name")
      set(BLAS_INCLUDE_SUFFIXES "extern/include")
      list(GET BLAS_LIBRARIES 0 _MBLAS_LIB)
      get_filename_component(_bdir ${_MBLAS_LIB} DIRECTORY)
      set(INCLUDE_DIR_HINTS "${_bdir}/../../")
    endif(BLAS_LIBRARIES)
  endif()

  ## Generic BLAS library? ##
  if((NOT BLAS_LIBRARIES)
      AND ((NOT WITH_BLAS) OR (WITH_BLAS STREQUAL "generic")))
    message(STATUS "Try to find a generic blas ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "blas"
      "")
    if(BLAS_LIBRARIES)
      set(WITH_BLAS "generic" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER cblas.h CACHE STRING "Blas header name")
    endif()
  endif()

  ## Generic BLAS library (if cblas not found in blas) ##
  if((NOT BLAS_LIBRARIES)
      AND ((NOT WITH_BLAS) OR (WITH_BLAS STREQUAL "generic")))
    message(STATUS "Try to find a generic blas (bis) ...")
    check_blas_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "cblas;blas"
      "")

    if(BLAS_LIBRARIES)
      set(WITH_BLAS "generic" CACHE STRING "Blas implementation type [mkl/openblas/atlas/accelerate/generic]" FORCE)
      set(BLAS_HEADER cblas.h CACHE STRING "Blas header name")
    endif()
  endif()

  ## If everything goes well, we have BLAS_LIBRARIES and some hints on where to find the headers.
  if(BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
  else(BLAS_LIBRARIES)
    set(BLAS_FOUND FALSE)
  endif(BLAS_LIBRARIES)

  ## Now the headers ...
  if(BLAS_FOUND)
    message(STATUS "Blas libraries have been found : ${BLAS_LIBRARIES}. We now turn to headers.")
    unset(BLAS_FOUND CACHE)
    # Get an hint for the location of blas headers : probably BLAS_LIBRARIES/../include or BLAS_DIR if given
    list(GET BLAS_LIBRARIES 0 BLAS_LIB)
    get_filename_component(_bdir ${BLAS_LIB} PATH)
    set(BLAS_LIBRARY_DIR ${_bdir} CACHE PATH "Blas libraries location." FORCE)

    # Either BLAS_DIR has been set by user or we use BLAS_LIBRARIES location to set it.
    if(NOT BLAS_DIR)
      get_filename_component(_bdir ${BLAS_LIBRARY_DIR} PATH)
      set(BLAS_DIR ${_bdir} CACHE PATH "Blas implementation location." FORCE)
    endif()

    set(BLAS_INC_DIR ${BLAS_DIR}/include)
    ## Use BLAS_DIR to set path for cmake search.
    ## If pkg-config has been used, maybe it gives some hints about headers location ...
    message(STATUS pkg hints : ${INCLUDE_DIR_HINTS})
    message(STATUS BLAS_INC_DIR : ${BLAS_INC_DIR})
    ## Note Franck : it seems that find_path process does not work as expected on Macosx : it does not respect the search sequence described
    # in cmake doc (i.e. CMAKE_PREFIX, then HINTS, PATHS ... ) and always turn to find apple framework cblas if
    # NO_DEFAULT_PATH is not set.
    if(APPLE) # First check in HINTS, no default, then global search.
      set(CMAKE_INCLUDE_PATH ${INCLUDE_DIR_HINTS} ${BLAS_INC_DIR})
      set(CMAKE_PREFIX_PATH ${BLAS_LIBRARIES})
      set(_headers)
      foreach(_file ${BLAS_HEADER})
	unset(_dir CACHE)
	find_path(_dir
	  NAMES ${_file}
	  HINTS ${BLAS_LIBRARIES}
	  PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
	  NO_DEFAULT_PATH
	  )
	find_path(_dir
	  NAMES ${_file}
	  PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
	  )
	if(_dir)
	  list(APPEND BLAS_INCLUDE_DIRS ${_dir})
	  list(APPEND _headers ${_dir}/${_file})
	 else()
	   if(BLAS_FIND_REQUIRED)
	     message(FATAL_ERROR "BLAS headers not found.")
	   endif()
	 endif()
      endforeach()
      unset(_dir CACHE)
      set(BLAS_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS} CACHE STRING "Blas headers location." FORCE)
    else() # The case which is supposed to always work
      find_path(BLAS_INCLUDE_DIRS
	NAMES ${BLAS_HEADER}
	PATH_SUFFIXES ${BLAS_INCLUDE_SUFFIXES}
        HINTS ${INCLUDE_DIR_HINTS}  ${BLAS_INC_DIR}
	)
    endif()
    if(BLAS_INCLUDE_DIRS)
      set(BLAS_FOUND 1 CACHE BOOL "Blas lib and headers have been found.")
    endif()
  endif()

  if(BLAS_FOUND) # SiconosConfig.h setup
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")

    if(WITH_BLAS STREQUAL "mkl")
      set(HAS_MKL_CBLAS 1 CACHE BOOL "Blas comes from Intel MKL.")

    elseif(WITH_BLAS STREQUAL "accelerate")
      set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")

    elseif(WITH_BLAS STREQUAL "atlas")
      set(HAS_ATLAS_CBLAS 1 CACHE BOOL "Blas  comes from Atlas framework ")

    elseif(WITH_BLAS STREQUAL "openblas")
      set(HAS_OpenBLAS 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")

    else()
      set(HAS_GenericCBLAS 1 CACHE BOOL "Blas is available from an unknown version.")

    endif()
  endif()
  set(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "Blas libraries.")

  if(NOT BLAS_FOUND AND BLAS_FIND_REQUIRED)
    message(FATAL_ERROR "Cannot find a library with BLAS API. Please specify library location using BLAS_DIR option or set your environment variables properly.")
  endif (NOT BLAS_FOUND AND BLAS_FIND_REQUIRED)
  if(NOT BLAS_FIND_QUIETLY)
    if(BLAS_FOUND)
      message(STATUS "Found a library with BLAS API (${WITH_BLAS}).")
      message("    Blas libraries : ${BLAS_LIBRARIES}.")
      message("    Blas headers : ${BLAS_HEADER} in ${BLAS_INCLUDE_DIRS}.")
    else(BLAS_FOUND)
      message(STATUS "Cannot find a library with BLAS API. Maybe you can try again using BLAS_DIR option or set your environment variables properly.")
    endif(BLAS_FOUND)
  endif(NOT BLAS_FIND_QUIETLY)

endif()
