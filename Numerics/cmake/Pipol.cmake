# On a Pipol system, set:
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

SET(PIPOL_USER $ENV{PIPOL_USER})
IF(NOT PIPOL_USER)
  SET(PIPOL_USER $ENV{USER})
  IF(NOT PIPOL_USER)
    MESSAGE(STATUS "neither $PIPOL_USER nor $USER exists")
  ENDIF(NOT PIPOL_USER)
ENDIF(NOT PIPOL_USER)


# scp is needed to get pipol-sub
FIND_PROGRAM(HAVE_SCP scp)
FIND_PROGRAM(HAVE_SSH ssh)


IF(PIPOL_USER)
  IF(HAVE_SSH)
    # get pipol systems
    EXECUTE_PROCESS(COMMAND 
      ssh ${PIPOL_USER}@pipol.inria.fr pipol-sub --query=systems 
      OUTPUT_VARIABLE PIPOL_SYSTEMS OUTPUT_STRIP_TRAILING_WHITESPACE)
  ENDIF(HAVE_SSH)

  IF(HAVE_SCP)
    MACRO(PIPOL_TARGET
        SYSTEM_PATTERN)

      # defaults
      IF(NOT PIPOL_CONFIGURE_COMMAND)
        SET(PIPOL_CONFIGURE_COMMAND "cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} ${CMAKE_SOURCE_DIR}")
      ENDIF(NOT PIPOL_CONFIGURE_COMMAND)

      IF(NOT PIPOL_MAKE_COMMAND)
        SET(PIPOL_MAKE_COMMAND make -j 2 -$(MAKEFLAGS))
      ENDIF(NOT PIPOL_MAKE_COMMAND)

      IF(NOT PIPOL_MAKE_TEST_COMMAND)
        SET(PIPOL_MAKE_TEST_COMMAND make -$(MAKEFLAGS) test)
      ENDIF(NOT PIPOL_MAKE_TEST_COMMAND)

      IF(NOT PIPOL_MAKE_INSTALL_COMMAND)
        SET(PIPOL_MAKE_INSTALL_COMMAND sudo make -$(MAKEFLAGS) install)
      ENDIF(NOT PIPOL_MAKE_INSTALL_COMMAND)

      ADD_CUSTOM_TARGET(
        ${SYSTEM_PATTERN}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND scp -q pipol.inria.fr:/usr/local/bin/pipol-sub . \;\\
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} ${PIPOL_RC_DIR_OPTION}  
        \"sudo mkdir -p /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          sudo chown \${USER} /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          cd /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          ${PIPOL_CONFIGURE_COMMAND} \;
          ${PIPOL_MAKE_COMMAND} \;
          ${PIPOL_MAKE_TEST_COMMAND} \; 
          ${PIPOL_MAKE_INSTALL_COMMAND} \"
        )

      ADD_CUSTOM_TARGET(
        package-${SYSTEM_PATTERN}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND scp -q pipol.inria.fr:/usr/local/bin/pipol-sub . \;\\
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} ${PIPOL_RC_DIR_OPTION}  
        \"sudo mkdir -p /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          sudo chown \${USER} /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          cd /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} \;
          ${PIPOL_CONFIGURE_COMMAND} \;
          ${PIPOL_MAKE_COMMAND} \;
          ${PIPOL_MAKE_TEST_COMMAND} \; 
          ${PIPOL_MAKE_INSTALL_COMMAND} \"
        )
    ENDMACRO(PIPOL_TARGET)
  ENDIF(HAVE_SCP)
  
ENDIF(PIPOL_USER)

IF(PIPOL_RC_DIR)
  SET(PIPOL_RC_DIR_OPTION --rc-dir=${PIPOL_RC_DIR})
ENDIF(PIPOL_RC_DIR)

# add a target for each pipol system
IF(PIPOL_SYSTEMS)
  FOREACH(SYSTEM ${PIPOL_SYSTEMS})
    # target with rc-dir
    PIPOL_TARGET(${SYSTEM})
  ENDFOREACH(SYSTEM ${PIPOL_SYSTEMS})
ENDIF(PIPOL_SYSTEMS)
