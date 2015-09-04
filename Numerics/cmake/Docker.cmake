# Compilations and tests inside docker containers

# DOCKER_IMAGE specifications in ${CMAKE_SOURCE_DIR}/../Build/Docker/${DOCKER_IMAGE}
# DOCKER_SHARED_DIRECTORIES : directories visibles in the container (source directory and binary directory are visibles)
# DOCKER_CMAKE_FLAGS : flags passed to cmake 
# DOCKER_MAKE_FLAGS  : flags passed to make
# DOCKER_MAKE_TEST_FLAGS  : flags passed to make test

macro(add_docker_targets)

  set(options)
  set(oneValueArgs DOCKER_USER DOCKER_UID DOCKER_GID DOCKER_TEMPLATE DOCKER_IMAGE DOCKER_IMAGE_DIR DOCKER_REPOSITORY DOCKER_FILE DOCKER_WORKDIR DOCKER_HOST_INSTALL_PREFIX DOCKER_SHARED_DIRECTORIES DOCKER_CMAKE_FLAGS DOCKER_MAKE_FLAGS DOCKER_MAKE_TEST_FLAGS)

  if (NOT DOCKER_USER)
    set(DOCKER_USER $ENV{USER})
  endif()
  
  if (NOT DOCKER_UID)
    execute_process(COMMAND id -u OUTPUT_VARIABLE _DOCKER_UID)
    string(STRIP ${_DOCKER_UID} DOCKER_UID)
  endif()
  
  if (NOT DOCKER_GID)
    execute_process(COMMAND id -g OUTPUT_VARIABLE _DOCKER_GID)
    string(STRIP ${_DOCKER_GID} DOCKER_GID)
  endif()
  
  if(DOCKER_TEMPLATE)
    set(DOCKER_IMAGE ${DOCKER_DISTRIB}-${DOCKER_TEMPLATE})
    set(DOCKER_FILE ${DOCKER_IMAGE_DIR}/${DOCKER_TEMPLATE})
  endif()
  
  if(NOT DOCKER_IMAGE)
    message(FATAL_ERROR "Docker : DOCKER_IMAGE unset")
  else()
    string(REPLACE ":" "-" DOCKER_IMAGE_AS_DIR ${DOCKER_IMAGE})
    string(REPLACE "." "-" DOCKER_IMAGE_AS_DIR ${DOCKER_IMAGE_AS_DIR})
  endif()
  
  if(NOT DOCKER_IMAGE_DIR)
    message(FATAL_ERROR "Docker : DOCKER_IMAGE_DIR unset")
  endif()
  
  if(NOT DOCKER_REPOSITORY)
    message(FATAL_ERROR "Docker : DOCKER_REPOSITORY unset")
  endif()
  
  if(NOT DOCKER_FILE)
    set(DOCKER_FILE ${DOCKER_IMAGE_DIR}/${DOCKER_IMAGE})
  endif()
  
  if(NOT DOCKER_WORKDIR)
    set(DOCKER_WORKDIR ${CMAKE_BINARY_DIR}/Docker/${DOCKER_IMAGE_AS_DIR})
  endif()
  file(MAKE_DIRECTORY ${DOCKER_WORKDIR})
  
  if(NOT DOCKER_HOST_INSTALL_PREFIX)
    set(DOCKER_HOST_INSTALL_PREFIX /tmp/$ENV{USER}/docker/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR})
    file(MAKE_DIRECTORY ${DOCKER_HOST_INSTALL_PREFIX})
  endif()
  
  message(STATUS "Docker cmake flags : ${DOCKER_CMAKE_FLAGS}")
  message(STATUS "Docker make flags : ${DOCKER_MAKE_FLAGS}")
  message(STATUS "Docker make test flags : ${DOCKER_MAKE_TEST_FLAGS}")
  message(STATUS "Docker make test flags : ${DOCKER_MAKE_TEST_FLAGS}")
  
  set(GENERATED_DOCKER_FILE ${CMAKE_CURRENT_BINARY_DIR}/Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR}/Dockerfile)
  configure_file(${DOCKER_FILE} ${GENERATED_DOCKER_FILE})
  
  file(APPEND ${GENERATED_DOCKER_FILE} "RUN groupadd --gid=${DOCKER_GID} ${DOCKER_USER}\n")
  file(APPEND ${GENERATED_DOCKER_FILE} "RUN useradd --gid=${DOCKER_GID} --uid=${DOCKER_UID} --create-home --shell=/bin/sh ${DOCKER_USER}\n")
  file(APPEND ${GENERATED_DOCKER_FILE} "USER ${DOCKER_USER}\n")
  file(APPEND ${GENERATED_DOCKER_FILE} "WORKDIR /home/${DOCKER_USER}\n")
  
  set(DOCKER_VFLAGS)
  foreach(_D ${DOCKER_SHARED_DIRECTORIES};${CMAKE_SOURCE_DIR};${CMAKE_BINARY_DIR})
    set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${_D}:${_D})
  endforeach()
  
  set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${DOCKER_HOST_INSTALL_PREFIX}:${CMAKE_INSTALL_PREFIX})
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-build
    COMMENT "Docker Build : ${DOCKER_IMAGE}"
    COMMAND cd Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR} && docker build -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} .
    )
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-cmake
    COMMENT "Docker cmake : ${DOCKER_IMAGE}"
    COMMAND docker run --user=${DOCKER_USER} --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} cmake ${CMAKE_SOURCE_DIR} ${DOCKER_CMAKE_FLAGS})
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make
    COMMENT "Docker make : ${DOCKER_IMAGE}"
    COMMAND docker run --user=${DOCKER_USER} --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_FLAGS})
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-test
    COMMENT "Docker make test : ${DOCKER_IMAGE}"
    COMMAND docker run --user=${DOCKER_USER} --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} test)
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-install
    COMMENT "Docker make install : ${DOCKER_IMAGE}"
    COMMAND docker run --user=${DOCKER_USER} --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_INSTALL_FLAGS} install)
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-clean
    COMMENT "Docker make clean : ${DOCKER_IMAGE}"
    COMMAND docker run --user=${DOCKER_USER} --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_CLEAN_FLAGS} clean)
  
endmacro()
