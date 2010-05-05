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
#    EXECUTE_PROCESS(COMMAND ssh ${PIPOL_USER}@pipol.inria.fr pipol-sub si OUTPUT_VARIABLE PIPOL_SYSTEMS OUTPUT_STRIP_TRAILING_WHITESPACE)
# WHY ??
SET(PIPOL_SYSTEMS amd64_kvm-linux-debian-lenny;amd64_kvm-windows-7;amd64-linux-centos-5;amd64-linux-debian-etch;amd64-linux-debian-lenny;amd64-linux-debian-testing;amd64-linux-fedora-core10;amd64-linux-fedora-core11;amd64-linux-fedora-core12;amd64-linux-fedora-core7;amd64-linux-fedora-core8;amd64-linux-fedora-core9;amd64-linux-mandriva-2007_springs_powerpack;amd64-linux-mandriva-2009_powerpack;amd64-linux-opensuse-11;amd64-linux-redhatEL-5.0;amd64-linux-suse-LES10;amd64-linux-ubuntu-feisty;amd64-linux-ubuntu-hardy;amd64-linux-ubuntu-intrepid;amd64-linux-ubuntu-jaunty;amd64-linux-ubuntu-karmic;amd64-linux-ubuntu-lucid;amd64-unix-freebsd-7;amd64-windows-server-2003-64bits;amd64-windows-server-2008-64bits;i386_kvm-linux-debian-lenny;i386_kvm-linux-fedora-core12;i386_kvm-windows-xp-pro-sp3;i386-linux-centos-5;i386-linux-debian-etch;i386-linux-debian-lenny;i386-linux-debian-testing;i386-linux-fedora-core10;i386-linux-fedora-core11;i386-linux-fedora-core12;i386-linux-fedora-core7;i386-linux-fedora-core8;i386-linux-fedora-core9;i386-linux-mandriva-2007_springs_powerpack;i386-linux-mandriva-2009_powerpack;i386-linux-opensuse-11;i386-linux-redhatEL-5.0;i386-linux-suse-LES10;i386-linux-ubuntu-feisty;i386-linux-ubuntu-hardy;i386-linux-ubuntu-intrepid;i386-linux-ubuntu-jaunty;i386-linux-ubuntu-karmic;i386-linux-ubuntu-lucid;i386_mac-mac-osx-server-leopard;i386-unix-freebsd-7;i386-unix-opensolaris-10;i386-unix-opensolaris-11;i386-unix-solaris-10;ia64-linux-debian-lenny.dd;ia64-linux-fedora-core9.dd;ia64-linux-redhatEL-5.0.dd;x86_64_mac-mac-osx-server-snow-leopard;x86_mac-mac-osx-server-snow-leopard)
  ENDIF(HAVE_SSH)

  IF(HAVE_SCP)
    MACRO(PIPOL_TARGET
        SYSTEM_PATTERN
        RC_DIR)

      IF(NOT PIPOL_CONFIGURE_COMMAND)
        SET(PIPOL_CONFIGURE_COMMAND cmake ${CMAKE_SOURCE_DIR})
      ENDIF(NOT PIPOL_CONFIGURE_COMMAND)

      IF(NOT PIPOL_MAKE_COMMAND)
        SET(PIPOL_MAKE_COMMAND make -j 2)
      ENDIF(NOT PIPOL_MAKE_COMMAND)

      IF(NOT PIPOL_MAKE_TEST_COMMAND)
        SET(PIPOL_MAKE_TEST_COMMAND make test)
      ENDIF(NOT PIPOL_MAKE_TEST_COMMAND)

      IF(NOT PIPOL_MAKE_INSTALL_COMMAND)
        SET(PIPOL_MAKE_INSTALL_COMMAND sudo make install)
      ENDIF(NOT PIPOL_MAKE_INSTALL_COMMAND)

      ADD_CUSTOM_TARGET(
        ${SYSTEM_PATTERN}
        COMMENT "PIPOL Build : ${SYSTEM_PATTERN}"
        COMMAND scp -q pipol.inria.fr:/usr/local/bin/pipol-sub .
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 --export=${CMAKE_SOURCE_DIR} --rc-dir=${RC_DIR} mkdir -p /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME}
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} ${PIPOL_CONFIGURE_COMMAND}
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} ${PIPOL_MAKE_COMMAND}
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} ${PIPOL_MAKE_TEST_COMMAND}
        COMMAND ./pipol-sub --pipol-user=${PIPOL_USER} ${SYSTEM_PATTERN} 02:00 --reconnect --group --keep --verbose=1 cmake -E chdir /pipol/\${USER}/${CMAKE_BUILD_TYPE}/${PROJECT_SHORT_NAME} ${PIPOL_MAKE_INSTALL_COMMAND}
        )
    ENDMACRO(PIPOL_TARGET)
  ENDIF(HAVE_SCP)
  
ENDIF(PIPOL_USER)
