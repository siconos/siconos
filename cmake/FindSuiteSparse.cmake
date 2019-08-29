#  find_package(SuiteSparse COMPONENTS CSparse)

# - Try to find SuiteSparse
# Once done this will define
#  SuiteSparse_FOUND - System has SuiteSparse
#  SuiteSparse_INCLUDE_DIRS - The SuiteSparse include directories
#  SuiteSparse_LIBRARIES - The libraries needed to use SuiteSparse

# For each component (currently only "CXSparse"):
#  SuiteSparse_(component)_FOUND
#  SuiteSparse_(component)_LIBRARY
#  SuiteSparse_(component)_INCLUDE_DIR
#  SuiteSparse_(component)_DEFINITIONS


# Required :
# - header : cs.h
# - libs : colamd, cxsparse (for cxsparse)

include(FindPackageHandleStandardArgs)


# Provide SuiteSparse_<C> variables for each component.
foreach(component IN LISTS SuiteSparse_FIND_COMPONENTS)
  set(SuiteSparse_USE_${component} 1 )
endforeach()

set(_SUITESPARSE_REQUIRED_VARS)

# -- cxsparse component --
if (SuiteSparse_USE_CXSparse)
  find_path(SuiteSparse_CXSparse_INCLUDE_DIR cs.h
    PATH_SUFFIXES SuiteSparse suitesparse
    DOC "Directory containing CXSparse header")

  find_library(SuiteSparse_CXSparse_LIBRARY NAMES cxsparse)
  # if (CXSparse_INCLUDE_DIR AND CXSparse_LIBRARY)
  #     set(SuiteSparse_CXSparse_FOUND TRUE)
  # endif()

  # At least on some systems we need to link to libcolamd which is
  # another output from suitesparse.
  find_library(colamd_LIBRARY NAMES colamd)

  list(APPEND _SUITESPARSE_REQUIRED_VARS
    SuiteSparse_CXSparse_LIBRARY SuiteSparse_CXSparse_INCLUDE_DIR)  
  set(SuiteSparse_CXSparse_LIBRARIES ${SuiteSparse_CXSparse_LIBRARY} ${colamd_LIBRARY})
  if(SuiteSparse_CXSparse_INCLUDE_DIR AND  SuiteSparse_CXSparse_LIBRARY)
    set(SuiteSparse_CXSparse_FOUND TRUE)
  endif()
endif ()

if(_SUITESPARSE_REQUIRED_VARS)
  find_package_handle_standard_args(SuiteSparse
    REQUIRED_VARS ${_SUITESPARSE_REQUIRED_VARS}
    HANDLE_COMPONENTS)
else()
  set(SuiteSparse_FOUND)
endif()

foreach(_component IN LISTS SuiteSparse_FIND_COMPONENTS)
  if(SuiteSparse_${_component}_FOUND AND NOT TARGET SuiteSparse::${_component})
    add_library(SuiteSparse::${_component} UNKNOWN IMPORTED)
    set_target_properties(SuiteSparse::${_component} PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION "${SuiteSparse_${_component}_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${SuiteSparse_${_component}_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${SuiteSparse_${_component}_LIBRARIES}")
    mark_as_advanced(SuiteSparse_${_component}_LIBRARIES SuiteSparse_${_component}_INCLUDE_DIR)
  endif()
endforeach()


