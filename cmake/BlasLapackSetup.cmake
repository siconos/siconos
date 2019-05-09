# - Try to find Blas (cblas) and a complete Lapack implementation.
# Once done, this will define
#
#  HAS_CBLAS - system has CBLAS
#  BLAS_INCLUDE_DIRS - The directories  which contain cblas headers.
#  BLAS_LIBRARIES - Link these to use cblas.
#
# Depending on the implementation of cblas, another var is set :
# - HAS_MKL_CBLAS for mkl intel
# - HAS_OpenBLAS for openblas (ex-goto)
# - HAS_ACCELERATE for apple accelerate framework
# - HAS_ATLAS_CBLAS for atlas cblas. 
# 
#  For Lapack, we set the following :
#  LAPACK_INCLUDE_DIRS  - The directories which contain lapack headers.
#  LAPACK_LIBRARIES - Link these to use lapack.
# Depending on the implementation of cblas, another var is set :
#  - HAS_MKL_LAPACKE for lapacke in mkl intel
#  - HAS_ATLAS_LAPACK for (uncomplete) lapack from atlas
#  - HAS_LAPACKE for any lapack distribution based on lapacke. See http://www.netlib.org/lapack/lapacke.html.
#
# User-defined options :
# - USE_MKL to force mkl intel
# - USE_ATLAS to force atlas
# - USE_APPLE_FRAMEWORK to force apple framework
# - USE_OpenBLAS to force openblas
# 
# To find a specific (local) blas and/or lapack, set LD_LIBRARY_PATH (or equivalent) properly. 
#
# 
# inspired from http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries
# 
# Note : the following implementations have been succesfully tested for Siconos :
#  - netlib blas/cblas + lapacke 
#  - mkl cblas + lapacke
#  - OpenBlas (GotoBlas "offspring")
#  - Atlas 
#


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

#==================
#    CBlas and Lapack
#==================
# We set BLA_VENDOR = Generic as a default value to avoid BLA_VENDOR=All
if(NOT BLA_VENDOR)
  set(BLA_VENDOR Generic)
endif()

set(WITH_CBLAS 1)

if(USE_MKL)
 if(NOT BLA_VENDOR)
  set(BLA_VENDOR Intel10_64lp_seq)
 endif()
  unset(WITH_CBLAS)
elseif(USE_ATLAS)
  set(BLA_VENDOR ATLAS)
  if(NOT LAPACK_VENDOR)
    set(LAPACK_VENDOR "Generic")
  endif()
elseif(USE_OpenBLAS)
  set(BLA_VENDOR OpenBLAS)
  if(NOT LAPACKE_VENDOR)
    set(LAPACKE_VENDOR "Generic")
  endif()
elseif(USE_APPLE_FRAMEWORK)
  set(BLA_VENDOR Apple)
elseif(USE_MATLAB)
  set(BLA_VENDOR MATLAB)
  # TODO: LAPACK or LAPACKE
endif()

## this finds both BLAS and LAPACK
find_package(LAPACK REQUIRED)

#==================
#    CBlas
#==================
list(GET BLAS_LIBRARIES 0 BLAS_LIB)
get_filename_component(BLAS_DIR ${BLAS_LIB} PATH)
set(BLAS_DIR ${BLAS_DIR}/..)

unset(BLAS_FOUND)
## --- mkl blas ---
if(BLAS_LIBRARIES MATCHES "mkl.*")
  message(STATUS "Check for blas headers in mkl ...")
  find_path(BLAS_INCLUDE_DIRS
    NAMES mkl_cblas.h
    PATHS ${BLAS_DIR}/..
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    )
  
  if(BLAS_INCLUDE_DIRS)
    set(BLAS_FOUND 1)
    set(HAS_MKL_CBLAS 1 CACHE BOOL "Blas comes from Intel MKL.")
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  endif()
## --- Apple framework blas ---
elseif(BLAS_LIBRARIES MATCHES "Accelerate.*")
  message(STATUS "Try to find apple headers for blas ...")
  find_path(CBLAS_INCLUDE_DIR
    NAMES cblas.h
    )
  find_path(ACC_INCLUDE_DIR
    NAMES Accelerate.h
    )
  set(BLAS_INCLUDE_DIRS ${CBLAS_INCLUDE_DIR} ${ACC_INCLUDE_DIR})
  print_var(BLAS_INCLUDE_DIRS)
  if(BLAS_INCLUDE_DIRS)
    set(BLAS_FOUND 1)
    set(HAS_ACCELERATE 1 CACHE BOOL "Blas/Lapack come from Accelerate framework ")
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  endif()

## --- Atlas blas ---
elseif(BLAS_LIBRARIES MATCHES "atlas.*")
  find_path(BLAS_INCLUDE_DIRS
    NAMES cblas.h
    PATHS ${BLAS_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    )
  if(BLAS_INCLUDE_DIRS)
    set(BLAS_FOUND 1)
    set(HAS_ATLAS_CBLAS 1 CACHE BOOL "Blas  comes from Atlas framework ")
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  endif()

## --- OpenBlas (ex-Goto) ---
elseif(BLAS_LIBRARIES MATCHES "openblas.*")
  find_path(BLAS_INCLUDE_DIRS
    NAMES cblas.h
    PATHS ${BLAS_DIR}
    PATH_SUFFIXES include/openblas include
    NO_DEFAULT_PATH
    )
  if(BLAS_INCLUDE_DIRS)
    set(BLAS_FOUND 1)
    set(HAS_OpenBLAS 1 CACHE BOOL "Blas/Lapack come from OpenBlas ")
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  endif()
## - default case ##
else()
 find_path(BLAS_INCLUDE_DIRS
    NAMES cblas.h
    PATHS ${BLAS_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
    )
 find_path(BLAS_INCLUDE_DIRS
    NAMES cblas.h
    )
  if(BLAS_INCLUDE_DIRS)
    set(BLAS_FOUND 1)
    set(HAS_GenericCBLAS 1 CACHE BOOL "Blas is available from an unknown version.")
    set(HAS_CBLAS 1 CACHE BOOL "A CBlas implementation is available.")
  endif()
  
endif()

include_directories(${BLAS_INCLUDE_DIRS})
remember_link_libraries("${BLAS_LIBRARIES}")

message(STATUS "Blas Libraries : ${BLAS_LIBRARIES}")
message(STATUS "Blas include : ${BLAS_INCLUDE_DIRS}")
message(STATUS "Blas linker flags : ${BLAS_LINKER_FLAGS}")
message(STATUS "Blas found : ${BLAS_FOUND}")
message(STATUS "Blas Vendor : ${BLA_VENDOR}")

#==================
#    Lapack
#==================


# Apple framework 
if(HAS_ACCELERATE)
  # Note that if blas commes from accelerate, the path to lapack headers
  #  is the same as the one for blas headers.
  find_path(LAPACK_INCLUDE_DIRS
    NAMES clapack.h
    PATHS ${BLAS_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    )
  set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/clapack.h)
  set(LAPACK_SUFFIX "_")
  set(LAPACK_PREFIX)

else()
  # Get location of blas header as hint to check for lapack headers.
  list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
  get_filename_component(LAPACK_DIR ${LAPACK_LIB} PATH)
  set(LAPACK_DIR ${LAPACK_DIR}/..)

  unset(LAPACK_FOUND)
  # mkl intel 
  if(HAS_MKL_CBLAS)
    message(STATUS "Try to find lapack headers in mkl ...")
    find_path(LAPACK_INCLUDE_DIRS
      NAMES mkl_lapacke.h
      PATHS ${LAPACK_DIR}/..
      PATH_SUFFIXES include
      NO_DEFAULT_PATH
      )
    
    if(LAPACK_INCLUDE_DIRS)
      set(LAPACK_FOUND 1)
      set(HAS_MKL_LAPACKE 1 CACHE BOOL "LAPACK comes from Intel MKL.")
      set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/mkl_lapacke.h)
      set(LAPACK_SUFFIX)
      set(LAPACK_PREFIX "LAPACKE_")
    endif()
  # we can have openblas and lapack from atlas
  elseif(HAS_ATLAS_CBLAS)
    message("Try to find lapack headers in atlas ...")
    find_path(LAPACK_INCLUDE_DIRS
      NAMES clapack.h
      PATHS ${LAPACK_DIR}
      PATH_SUFFIXES include include/atlas
      NO_DEFAULT_PATH
      )
    if(LAPACK_INCLUDE_DIRS)
      set(LAPACK_FOUND 1)
      set(HAS_ATLAS_LAPACK 1 CACHE BOOL "LAPACK comes from atlas.")
      set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/clapack.h)
      set(LAPACK_SUFFIX)
      set(LAPACK_PREFIX "clapack_")
    endif()
  else() # netlib, openblas
    find_path(LAPACK_INCLUDE_DIRS
      NAMES lapacke.h
      PATHS ${LAPACK_DIR}
      PATH_SUFFIXES include
      NO_DEFAULT_PATH
      )
    if(LAPACK_INCLUDE_DIRS)
      set(LAPACK_FOUND 1)
      set(HAS_LAPACKE 1 PARENT_SCOPE)
      set(LAPACK_HEADER ${LAPACK_INCLUDE_DIRS}/lapacke.h)
      set(LAPACK_SUFFIX)
      set(LAPACK_PREFIX "LAPACKE_")
    endif()
  endif()
endif()

include_directories("${LAPACK_INCLUDE_DIRS}")
remember_link_libraries("${LAPACK_LIBRARIES}")

message(STATUS "Lapack Libraries : ${LAPACK_LIBRARIES}")
message(STATUS "Lapack include : ${LAPACK_INCLUDE_DIRS}")
message(STATUS "Lapack linker flags : ${LAPACK_LINKER_FLAGS}")
message(STATUS "Lapack found : ${LAPACK_FOUND}")


# === Check for lapack functions ===
# We check only the functions that are known to be un-implemented
# in some lapack libraries (namely atlas ...)
# This is probably a temporary check since it's likely
# we will stop atlas checking for lapack? 

set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set(CMAKE_REQUIRED_INCLUDES ${BLAS_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
## dgesvd ##
unset(HAS_LAPACK_DGESVD)
set(GENERIC_NAME "DGESVD")
set(FUNC_NAME "${LAPACK_PREFIX}dgesvd${LAPACK_SUFFIX}")
check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGESVD)

## dgels ##
unset(HAS_LAPACK_DGELS)
set(GENERIC_NAME "DGELS")
set(FUNC_NAME "${LAPACK_PREFIX}dgels${LAPACK_SUFFIX}")
check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DGELS)

## dtrtrs ##
unset(HAS_LAPACK_DTRTRS)
set(GENERIC_NAME "DTRTRS")
set(FUNC_NAME "${LAPACK_PREFIX}dtrtrs${LAPACK_SUFFIX}")
check_lapack_has_function(${GENERIC_NAME} ${FUNC_NAME} ${LAPACK_HEADER} HAS_LAPACK_DTRTRS)



