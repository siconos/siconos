# ==========================================
# cmake utilities for gams library.
#
# See http://www.gams.com
# 
# Usage:
# include(gams_setup)
#==========================================

IF(GAMSCAPI_FOUND)

  LIST(APPEND ${COMPONENT}_SRCS ${GAMS_C_API_FILES})
  INCLUDE_DIRECTORIES(${GAMS_C_API_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES("${GAMS_DIR}/testlib_ml") # XXX Hack -- xhub
  INSTALL(DIRECTORY ${GAMS_MODELS_SOURCE_DIR}
   DESTINATION share/${PROJECT_NAME})

  INCLUDE(CheckCCompilerFlag)
  CHECK_C_COMPILER_FLAG("-Werror=conversion" C_HAVE_WERR_CONV)
  IF(C_HAVE_WERR_CONV)
   SET_SOURCE_FILES_PROPERTIES(${GAMS_C_API_FILES} PROPERTIES COMPILE_FLAGS
    "-Wno-error=conversion")
  ENDIF(C_HAVE_WERR_CONV)

  # XXX hack for now ...
  COMPILE_WITH(PathVI SICONOS_COMPONENTS numerics)
  IF(PathVI_FOUND)
    SET(HAVE_GAMS_PATHVI TRUE)
    #    compile_with(Cplex)
  ELSE(PathVI_FOUND)
    SET(HAVE_GAMS_PATHVI FALSE)
  ENDIF(PathVI_FOUND)

  # GAMS distributes Cplex ...
  #  IF(NOT Cplex_LIBRARY_DIRECTORY)
  #    SET(Cplex_LIBRARY_DIRECTORY ${GAMS_DIR})
  #   ENDIF()
 ENDIF(GAMSCAPI_FOUND)
