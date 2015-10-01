# ==========================================
# cmake utilities for gams library.
#
# See http://www.gams.com
# 
# Usage:
# include(gams_setup)
#==========================================

IF(GAMS_DIR)
  
  SET(GAMS_C_API_FIND_REQUIRED TRUE)
  COMPILE_WITH(GamsCApi)
  LIST(APPEND ${COMPONENT}_SRCS ${GAMS_C_API_FILES})
  INCLUDE_DIRECTORIES(${GAMS_C_API_INCLUDE_DIRS})
  INCLUDE_DIRECTORIES("${GAMS_DIR}/testlib_ml") # XXX Hack -- xhub
  SET(GAMS_MODELS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/share/gams")
  SET(GAMS_MODELS_SHARE_DIR "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/gams")
  INSTALL(DIRECTORY ${GAMS_MODELS_SOURCE_DIR}
   DESTINATION share/${PROJECT_NAME})

  INCLUDE(CheckCCompilerFlag)
  CHECK_C_COMPILER_FLAG("-Werror=conversion" C_HAVE_WERR_CONV)
  IF(C_HAVE_WERR_CONV)

   SET_SOURCE_FILES_PROPERTIES(${GAMS_C_API_FILES} PROPERTIES COMPILE_FLAGS
    "-Wno-error=conversion")
  ENDIF()

  # XXX hack for now ...
  COMPILE_WITH(PathVI)
  IF(HAVE_PATHVI)
    SET(HAVE_GAMS_PATHVI TRUE)
  ELSE(HAVE_PATHVI)
    SET(HAVE_GAMS_PATHVI FALSE)
  ENDIF(HAVE_PATHVI)

  # GAMS distributes Cplex ...
  IF(GAMS_FOUND AND NOT Cplex_LIBRARY_DIRECTORY)
    SET(Cplex_LIBRARY_DIRECTORY ${GAMS_DIR})
   ENDIF()
ENDIF(GAMS_DIR)
