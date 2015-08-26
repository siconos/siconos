# Compilations and tests inside docker containers

# DOCKER_IMAGE specifications in ${CMAKE_SOURCE_DIR}/../Build/Docker/${DOCKER_IMAGE}
# DOCKER_SHARED_DIRECTORIES : directories visibles in the container (source directory and binary directory are visibles)
# DOCKER_CMAKE_FLAGS : flags passed to cmake 
# DOCKER_MAKE_FLAGS  : flags passed to make
# DOCKER_MAKE_TEST_FLAGS  : flags passed to make test

if(NOT DOCKER_IMAGE)
  message(FATAL_ERROR "Docker : DOCKER_IMAGE unset")
endif()

if(NOT DOCKER_IMAGE_DIR)
  message(FATAL_ERROR "Docker : DOCKER_IMAGE_DIR unset")
endif()

if(NOT DOCKER_FILE)
  set(DOCKER_FILE ${DOCKER_IMAGE_DIR}/${DOCKER_IMAGE})
endif()

if(NOT DOCKER_WORKDIR)
  set(DOCKER_WORKDIR ${CMAKE_BINARY_DIR}/Docker/${DOCKER_IMAGE})
endif()

configure_file(${DOCKER_FILE} Docker/Context/Dockerfile)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/Docker/${DOCKER_IMAGE})

set(DOCKER_VFLAGS)
foreach(_D ${DOCKER_SHARED_DIRECTORIES};${CMAKE_SOURCE_DIR};${CMAKE_BINARY_DIR})
  set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${_D}:${_D})
endforeach()

set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${DOCKER_HOST_INSTALL_PREFIX}:${CMAKE_INSTALL_PREFIX})

ADD_CUSTOM_TARGET(
  docker-image
  COMMENT "Docker Build : ${DOCKER_IMAGE}"
  COMMAND cd Docker/Context && docker build -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} .
  )

ADD_CUSTOM_TARGET(
  docker-cmake
  COMMENT "Docker cmake : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} cmake ${CMAKE_SOURCE_DIR} ${DOCKER_CMAKE_FLAGS})

ADD_CUSTOM_TARGET(
  docker-make
  COMMENT "Docker make : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_FLAGS})

ADD_CUSTOM_TARGET(
  docker-make-test
  COMMENT "Docker make test : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} test)

ADD_CUSTOM_TARGET(
  docker-make-install
  COMMENT "Docker make test : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} install)

ADD_CUSTOM_TARGET(
  docker-make-clean
  COMMENT "Docker make clean : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true ${DOCKER_VFLAGS} --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} clean)
