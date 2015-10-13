# Find the Path Ferris includes and libraries.
# The following variables are set if Ferris is found.  If Ferris is not
# found, GAMS_FOUND is set to false.
#  GAMS_FOUND        - True when the GAMS include directory is found.
#  GAMS_C_API_INCLUDE_DIRS - the path to where the Path Ferris include files are.

INCLUDE(FindPackageHandleStandardArgs)

SET(FILES_TO_CHECK "idxcc.c;optcc.c;gamsxcc.c;gmomcc.c")

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
  SET(HAVE_GAMS_C_API TRUE BOOL "The GAMS C API has been found")
  SET(GAMS_FOUND TRUE CACHE BOOL "A GAMS install has been found")
ELSE(GAMS_C_API_OK)
  IF(GAMS_C_API_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR
      "Required GAMS C API not found. Please check the GAMS_DIR argument")
  ENDIF(GAMS_C_API_FIND_REQUIRED)
ENDIF(GAMS_C_API_OK)
