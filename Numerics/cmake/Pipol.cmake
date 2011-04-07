# Remote compilations and tests on Pipol
# PIPOL_USER must be set in the environment
#

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
# default
#

SET(PIPOL_USER $ENV{PIPOL_USER} CACHE STRING "user login on a pipol system")

SET(PIPOL_CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${CMAKE_SOURCE_DIR} CACHE STRING "the configuration command")

SET(PIPOL_CONFIGURE_COMMAND cmake -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${CMAKE_SOURCE_DIR} CACHE STRING "configuration command on a pipol system")

SET(PIPOL_MAKE_COMMAND_TARGET "" CACHE STRING "make command target on a pipol system")

SET(PIPOL_MAKE_COMMAND make -j 2 -$(MAKEFLAGS) CACHE STRING "make command on a pipol system")

SET(PIPOL_MAKE_TEST_COMMAND make -j 2 -$(MAKEFLAGS) test CACHE STRING "make test command on a pipol system")

SET(PIPOL_MAKE_INSTALL_COMMAND sudo make -$(MAKEFLAGS) install CACHE STRING "make install command on a pipol system")

SET(PIPOL_POST_INSTALL_COMMAND sudo chown ${PIPOL_USER} install_manifest.txt  CACHE STRING "post install command on a pipol system")

SET(PIPOL_PACKAGE_COMMAND cpack -G \\\$$PIPOL_CPACK_G . CACHE STRING "package command on a pipol system")

SET(PIPOL_SUB_RSYNC_OPTIONS -aC CACHE STRING "rsync options")

SET(PIPOL_DURATION 02:00 CACHE STRING "duration of system reservation")

# ssh/rsync mandatory 
FIND_PROGRAM(HAVE_SSH ssh)
FIND_PROGRAM(HAVE_RSYNC rsync)

IF(PIPOL_USER)
  MESSAGE(STATUS "Pipol user is ${PIPOL_USER}")
  IF(HAVE_SSH)
    # get pipol systems
    EXECUTE_PROCESS(COMMAND 
      ssh ${PIPOL_USER}@pipol.inria.fr pipol-sub --query=systems 
      OUTPUT_VARIABLE PIPOL_SYSTEMS OUTPUT_STRIP_TRAILING_WHITESPACE)
  ENDIF(HAVE_SSH)

  IF(HAVE_RSYNC)
    MACRO(PIPOL_TARGET
        SYSTEM_PATTERN)

      STRING(REPLACE ".dd.gz" "" SYSTEM_TARGET ${SYSTEM_PATTERN})

      ADD_CUSTOM_TARGET(
        ${SYSTEM_TARGET}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND rsync ${PIPOL_USER}@pipol.inria.fr:/usr/local/bin/pipol-sub . 
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} ${PIPOL_DURATION} --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} ${PIPOL_RC_DIR_OPTION} --rsynco=${PIPOL_SUB_RSYNC_OPTIONS}
        \"sudo mkdir -p \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          sudo chown ${PIPOL_USER} \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          cd \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          ${PIPOL_CONFIGURE_COMMAND} \;
          ${PIPOL_MAKE_COMMAND} ${PIPOL_MAKE_COMMAND_TARGET}\;
          ${PIPOL_MAKE_TEST_COMMAND} \; 
          ${PIPOL_MAKE_INSTALL_COMMAND} \;
          ${PIPOL_POST_INSTALL_COMMAND} \"
        )

      ADD_CUSTOM_TARGET(
        make-${SYSTEM_TARGET}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND rsync ${PIPOL_USER}@pipol.inria.fr:/usr/local/bin/pipol-sub . 
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} ${PIPOL_DURATION} --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} ${PIPOL_RC_DIR_OPTION} --rsynco=${PIPOL_SUB_RSYNC_OPTIONS}
        \"sudo mkdir -p \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          sudo chown ${PIPOL_USER} \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          cd \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          ${PIPOL_MAKE_COMMAND} ${PIPOL_MAKE_COMMAND_TARGET}\"
        )

      ADD_CUSTOM_TARGET(
        package-${SYSTEM_TARGET}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND rsync ${PIPOL_USER}@pipol.inria.fr:/usr/local/bin/pipol-sub . \;\\
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} ${PIPOL_DURATION} --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} ${PIPOL_RC_DIR_OPTION}  
        \"sudo mkdir -p \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
        sudo chown ${PIPOL_USER} \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
        cd \\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME} \;
          ${PIPOL_CONFIGURE_COMMAND} \;
          ${PIPOL_MAKE_COMMAND} ${PIPOL_MAKE_COMMAND_TARGET}\;
          ${PIPOL_PACKAGE_COMMAND} \" \;\\
       COMMAND rsync -av ${PIPOL_USER}@`./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} ${PIPOL_DURATION} --reconnect --group --keep --verbose=0 echo \\\\$$PIPOL_HOST`.inrialpes.fr:`./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} ${PIPOL_DURATION} --reconnect --group --keep --verbose=0 ls \\\\$$PIPOL_WDIR/${PIPOL_USER}/${CMAKE_BUILD_TYPE}/${PROJECT_NAME}/*.\\\\$$PIPOL_PACK_EXT` .
       )
    ENDMACRO(PIPOL_TARGET)
  ENDIF(HAVE_RSYNC)
  
ENDIF(PIPOL_USER)

IF(PIPOL_RC_DIR)
  SET(PIPOL_RC_DIR_OPTION --rc-dir=${PIPOL_RC_DIR})
ENDIF(PIPOL_RC_DIR)

# add a target for each pipol system
IF(PIPOL_SYSTEMS)
  MESSAGE(STATUS "Adding Pipol targets")
  FOREACH(SYSTEM ${PIPOL_SYSTEMS})
    # target with rc-dir
    PIPOL_TARGET(${SYSTEM})
  ENDFOREACH(SYSTEM ${PIPOL_SYSTEMS})
ENDIF(PIPOL_SYSTEMS)
