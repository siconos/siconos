#  find_package(SuiteSparse COMPONENTS CSparse)

# - Try to find SuiteSparse
# Once done this will define
#  SuiteSparse_FOUND - System has SuiteSparse
#  SuiteSparse_INCLUDE_DIRS - The SuiteSparse include directories
#  SuiteSparse_LIBRARIES - The libraries needed to use SuiteSparse
#  SuiteSparse_DEFINITIONS - Compiler switches required for using SuiteSparse

# For each component (currently only "CXSparse"):
#  SuiteSparse_(component)_FOUND
#  SuiteSparse_(component)_LIBRARY
#  SuiteSparse_(component)_INCLUDE_DIR
#  SuiteSparse_(component)_DEFINITIONS

set(SuiteSparse_DEFINITIONS)

if( SuiteSparse_FIND_COMPONENTS )
  foreach( component ${SuiteSparse_FIND_COMPONENTS} )
    set( SuiteSparse_USE_${component} 1 )
  endforeach()
endif()

if (SuiteSparse_USE_CXSparse)
  find_path(CXSparse_INCLUDE_DIR cs.h
    PATH_SUFFIXES SuiteSparse suitesparse
    DOC "Directory containing CXSparse header")

  find_library(CXSparse_LIBRARY NAMES cxsparse)
  set(SuiteSparse_CXSparse_DEFINITIONS)
  if (CXSparse_INCLUDE_DIR AND CXSparse_LIBRARY)
      set(SuiteSparse_CXSparse_FOUND TRUE)
  endif()
endif ()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SuiteSparse_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SuiteSparse DEFAULT_MSG
                                  CXSparse_LIBRARY CXSparse_INCLUDE_DIR)

if (SuiteSparse_USE_CXSparse AND CXSparse_LIBRARY)
  mark_as_advanced(CXSparse_INCLUDE_DIR CXSparse_LIBRARY)
endif ()

set(SuiteSparse_LIBRARIES ${CXSparse_LIBRARY})
set(SuiteSparse_INCLUDE_DIRS ${CXSparse_INCLUDE_DIR})

# If at least one component was found, SuiteSparse was found.
# Users must check each component FOUND individually.
if( SuiteSparse_FIND_COMPONENTS )
  foreach( component ${SuiteSparse_FIND_COMPONENTS} )
    if (${component}_LIBRARY)
      set(SuiteSparse_FOUND TRUE)
    endif()
  endforeach()
endif()
