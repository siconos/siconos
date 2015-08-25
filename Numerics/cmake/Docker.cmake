# Compilations and tests inside docker containers

# DOCKER_IMAGE specifications in ${CMAKE_SOURCE_DIR}/../Build/Docker/${DOCKER_IMAGE}
# DOCKER_CMAKE_FLAGS : flags passed to cmake 
# DOCKER_MAKE_FLAGS  : flags passed to make
# DOCKER_MAKE_TEST_FLAGS  : flags passed to make test

if(NOT DOCKER_IMAGE)
  set(DOCKER_IMAGE debian-c-fortran-atlas)
endif()

set(DOCKER_FILE ${CMAKE_SOURCE_DIR}/../Build/Docker/${DOCKER_IMAGE})

configure_file(${DOCKER_FILE} DockerContext/Dockerfile)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/DockerBuild)

ADD_CUSTOM_TARGET(
  docker-image
  COMMENT "Docker Build : ${DOCKER_IMAGE}"
  COMMAND cd DockerContext && docker build -t siconos/${DOCKER_IMAGE} .
  )

ADD_CUSTOM_TARGET(
  docker-cmake
  COMMENT "Docker cmake : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true -v ${CMAKE_SOURCE_DIR}:${CMAKE_SOURCE_DIR} -v ${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR} --workdir=${CMAKE_BINARY_DIR}/DockerBuild -t siconos/${DOCKER_IMAGE} cmake ${CMAKE_SOURCE_DIR} ${DOCKER_CMAKE_FLAGS})

ADD_CUSTOM_TARGET(
  docker-make
  COMMENT "Docker make : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true -v ${CMAKE_SOURCE_DIR}:${CMAKE_SOURCE_DIR} -v ${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR} --workdir=${CMAKE_BINARY_DIR}/DockerBuild -t siconos/${DOCKER_IMAGE} make ${DOCKER_MAKE_FLAGS})

ADD_CUSTOM_TARGET(
  docker-make-test
  COMMENT "Docker make test : ${DOCKER_IMAGE}"
  COMMAND docker run --rm=true -v ${CMAKE_SOURCE_DIR}:${CMAKE_SOURCE_DIR} -v ${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR} --workdir=${CMAKE_BINARY_DIR}/DockerBuild -t siconos/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} test)
