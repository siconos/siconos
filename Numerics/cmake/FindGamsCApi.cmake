# Find the Path Ferris includes and libraries.
# The following variables are set if Ferris is found.  If Ferris is not
# found, PathFerris_FOUND is set to false.
#  PathFerris_FOUND        - True when the PathFerris include directory is found.
#  PathFerris_INCLUDE_DIRS - the path to where the Path Ferris include files are.
#  PathFerris_LIBRARY_DIRS - The path to where the Path library files are.
#  PathFerris_LIBRARIES    - The libraries to link against Path Ferris

# One may want to use a specific Ferris Library by setting
# PathFerris_LIBRARY_DIRECTORY before FIND_PACKAGE(PathFerris)
INCLUDE(FindPackageHandleStandardArgs)

SET(FILES_TO_CHECK "idxcc.c;optcc.c;gamsxcc.c;gevmcc.c")

SET(GAMS_C_API_DIR "apifiles/C/api")

SET(GAMS_C_API_OK TRUE)
SET(GAMS_C_API_FILES)

FOREACH(_F ${FILES_TO_CHECK})
  SET(_FF ${GAMS_DIR}/${GAMS_C_API_DIR}/${_F})
  IF(EXISTS ${_FF})
    LIST(APPEND GAMS_C_API_FILES ${_FF})
  ELSE(EXISTS ${_FF})
    MESSAGE(STATUS "The GAMS C API file ${_F} was not found (searched in ${GAMS_DIR}/${GAMS_C_API_DIR})")
    SET(GAMS_C_API_OK FALSE)
  ENDIF(EXISTS ${_FF})
ENDFOREACH(_F ${FILES_TO_CHECK})

IF(GAMS_C_API_OK)
  SET(GAMS_C_API_INCLUDE_DIRS ${GAMS_DIR}/${GAMS_C_API_DIR})
  SET(HAVE_GAMS_C_API 1)
ELSE(GAMS_C_API_OK)
  IF(GAMS_C_API_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required GAMS C API not found. Please check the GAMS_DIR argument")
  ENDIF(GAMS_C_API_FIND_REQUIRED)
ENDIF(GAMS_C_API_OK)
