# On Pipol system, set:
# PIPOL_IMAGE : Pipol image full name
# PIPOL_IMAGE_NAME : Pipol image without extension
# PIPOL_SITE : suggestion for SITE variable

SET(PIPOL_IMAGE $ENV{PIPOL_IMAGE})  
SET(_STMP "[${PIPOL_IMAGE}]")
IF(_STMP STREQUAL "[]")
  SET(PIPOL_IMAGE)
ELSE(_STMP STREQUAL "[]")
  GET_FILENAME_COMPONENT(PIPOL_IMAGE_NAME ${PIPOL_IMAGE} NAME_WE)
  SET(PIPOL_SITE "PIPOL")
ENDIF(_STMP STREQUAL "[]")

# 
#
#

MACRO(PIPOL_TARGET
    SYSTEM_PATTERN
    RC_DIR)

  ADD_CUSTOM_TARGET(
    ${SYSTEM_PATTERN}
    COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
    COMMAND scp -q pipol.inria.fr:/usr/local/bin/pipol-sub .
    COMMAND ./pipol-sub ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} --rc-dir=${RC_DIR} mkdir -p /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME}
  COMMAND ./pipol-sub ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} cmake ${CMAKE_SOURCE_DIR} -DWITH_TESTING=True
  COMMAND ./pipol-sub ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} make -j 2
  COMMAND ./pipol-sub ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} make test
  COMMAND ./pipol-sub ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} sudo make install
  )
ENDMACRO(PIPOL_TARGET)
