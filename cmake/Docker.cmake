# Compilations and tests inside docker containers

# DOCKER_IMAGE specifications in ${CMAKE_SOURCE_DIR}/Build/Docker/${DOCKER_IMAGE}
# DOCKER_SHARED_DIRECTORIES : directories visibles in the container (source directory and binary directory are visibles)
# DOCKER_CMAKE_FLAGS : flags passed to cmake 
# DOCKER_MAKE_FLAGS  : flags passed to make
# DOCKER_MAKE_TEST_FLAGS  : flags passed to make test

macro(add_docker_targets)

  set(options)
  set(oneValueArgs DOCKER_TEMPLATES DOCKER_TEMPLATE DOCKER_IMAGE DOCKER_IMAGE_DIR DOCKER_REPOSITORY DOCKER_FILE DOCKER_WORKDIR DOCKER_HOST_INSTALL_PREFIX DOCKER_SHARED_DIRECTORIES DOCKER_CMAKE_FLAGS DOCKER_MAKE_FLAGS DOCKER_MAKE_TEST_FLAGS DOCKER_CTEST_DRIVER)

  if(DOCKER_TEMPLATES)
    string(REPLACE "," ";" DOCKER_TEMPLATES_LIST ${DOCKER_TEMPLATES})

    string(REPLACE "," "-" DOCKER_IMAGE ${DOCKER_TEMPLATES})
    string(REPLACE "+" "x" DOCKER_IMAGE ${DOCKER_IMAGE})
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
    string(REPLACE "+" "x" DOCKER_IMAGE_AS_DIR ${DOCKER_IMAGE_AS_DIR})
  endif()
  
  if(NOT DOCKER_IMAGE_DIR)
    message(FATAL_ERROR "Docker : DOCKER_IMAGE_DIR unset")
  endif()
  
  if(NOT DOCKER_REPOSITORY)
    message(FATAL_ERROR "Docker : DOCKER_REPOSITORY unset")
  endif()
  
  set(GENERATED_DOCKER_FILE ${CMAKE_CURRENT_BINARY_DIR}/Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR}/Dockerfile)

  file(REMOVE ${GENERATED_DOCKER_FILE})
  
  message(STATUS "Docker templates list is : ${DOCKER_TEMPLATES_LIST}")

  foreach(_dt ${DOCKER_TEMPLATES_LIST})
    if(EXISTS ${DOCKER_IMAGE_DIR}/${_dt})
      file(READ ${DOCKER_IMAGE_DIR}/${_dt} _contents_)
      file(APPEND ${GENERATED_DOCKER_FILE} ${_contents_})
      list(REMOVE_ITEM DOCKER_TEMPLATES_LIST ${_dt})
    endif()
  endforeach()

  if(NOT DOCKER_MKSENV_SPLIT)
    set(DOCKER_MKSENV_SPLIT true)
  endif()

  execute_process(COMMAND ${DOCKER_MKSENV_COMMAND} --docker --distrib ${DOCKER_DISTRIB} --pkgs ${DOCKER_TEMPLATES} --split=${DOCKER_MKSENV_SPLIT} ${DOCKER_MKSENV_INPUT} OUTPUT_VARIABLE _contents_)
  file(APPEND ${GENERATED_DOCKER_FILE} ${_contents_})
  file(APPEND ${GENERATED_DOCKER_FILE} "RUN mkdir -p /usr/local\n")

  if(NOT DOCKER_WORKDIR)
    set(DOCKER_WORKDIR ${CMAKE_BINARY_DIR}/Docker/${DOCKER_IMAGE_AS_DIR})
  endif()
  file(MAKE_DIRECTORY ${DOCKER_WORKDIR})
  
  string(REPLACE "/" "-" DOCKER_WORKDIR_AS_NAME "workdir-${DOCKER_WORKDIR}")

  set(DOCKER_CMAKE_FLAGS_WITHOUT_DOCKER)
  foreach(_f ${DOCKER_CMAKE_FLAGS})
    string(REGEX MATCH "^-DDOCKER_.*" _mf ${_f})
    if (NOT _mf)
      list(APPEND DOCKER_CMAKE_FLAGS_WITHOUT_DOCKER ${_f})
    endif()
  endforeach()

  message(STATUS "Docker cmake flags : ${DOCKER_CMAKE_FLAGS_WITHOUT_DOCKER}")
  message(STATUS "Docker make flags : ${DOCKER_MAKE_FLAGS}")
  message(STATUS "Docker make test flags : ${DOCKER_MAKE_TEST_FLAGS}")

  set(DOCKER_VFLAGS)
  foreach(_D ${DOCKER_SHARED_DIRECTORIES};${CMAKE_SOURCE_DIR})
    set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${_D}:${_D})
  endforeach()
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-clean
    COMMENT "Docker clean : ${DOCKER_IMAGE}"
    COMMAND docker rm -f ${DOCKER_WORKDIR_AS_NAME}
    COMMAND docker rm -f ${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local
    COMMAND docker rmi -f ${DOCKER_REPOSITORY}/${DOCKER_IMAGE}
    )

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-build
    COMMENT "Docker Build : ${DOCKER_IMAGE}"
    COMMAND cd Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR} && docker build -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} .
    )

  add_custom_command(
    TARGET ${DOCKER_IMAGE_AS_DIR}-build
    PRE_BUILD 
    COMMENT "docker create workdir"
    COMMAND docker create --name=${DOCKER_WORKDIR_AS_NAME} -v ${DOCKER_WORKDIR} -i -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} true 2>/dev/null  && echo done || echo already done)
  
  add_custom_command(
    TARGET ${DOCKER_IMAGE_AS_DIR}-build
    PRE_BUILD 
    COMMENT "docker create /usr/local"
    COMMAND docker create --name=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local -v /usr/local -i -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} true 2>/dev/null && echo done || echo already done
    )
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-cmake
    COMMENT "Docker cmake : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} cmake ${CMAKE_SOURCE_DIR} ${DOCKER_CMAKE_FLAGS_WITHOUT_DOCKER})
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make
    COMMENT "Docker make : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_FLAGS})
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-test
    COMMENT "Docker make test : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} test)
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-install
    COMMENT "Docker make install : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_INSTALL_FLAGS} install)
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-clean
    COMMENT "Docker make clean : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_CLEAN_FLAGS} clean)

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-ctest
    COMMENT "Docker ctest : ${DOCKER_IMAGE}"
    COMMAND docker run --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_AS_NAME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} ctest -DCTEST_SOURCE_DIRECTORY=${CMAKE_SOURCE_DIR} -DCTEST_BINARY_DIRECTORY=${CMAKE_BINARY_DIR} -S ${DOCKER_CTEST_DRIVER} ${DOCKER_CMAKE_FLAGS_WITHOUT_DOCKER})

  if(NOT TARGET docker-clean)
    add_custom_target(
      docker-clean
      COMMENT "Docker clean"
      )
  endif()

  if(NOT TARGET docker-build)
    add_custom_target(
      docker-build ALL
      COMMENT "Docker build"
      )
  endif()

  if(NOT TARGET docker-cmake)
    add_custom_target(
      docker-cmake ALL
      COMMENT "Docker cmake"
      )
  endif()

  if(NOT TARGET docker-make)
    add_custom_target(
      docker-make ALL
      COMMENT "Docker make"
      )
  endif()

  if(NOT TARGET docker-make-test)
    add_custom_target(
      docker-make-test
      COMMENT "Docker make test"
      )
  endif()

  if(NOT TARGET docker-make-install)
    add_custom_target(
      docker-make-install
      COMMENT "Docker make install"
      )
  endif()

  if(NOT TARGET docker-make-clean)
    add_custom_target(
      docker-make-clean
      COMMENT "Docker make clean"
      )
  endif()

  if(NOT TARGET docker-ctest)
    add_custom_target(
      docker-ctest
      COMMENT "Docker ctest"
      )
  endif()

  add_dependencies(docker-clean ${DOCKER_IMAGE_AS_DIR}-clean)
  add_dependencies(docker-build ${DOCKER_IMAGE_AS_DIR}-build)
  add_dependencies(docker-cmake ${DOCKER_IMAGE_AS_DIR}-cmake)
  add_dependencies(docker-make ${DOCKER_IMAGE_AS_DIR}-make)
  add_dependencies(docker-make-test ${DOCKER_IMAGE_AS_DIR}-make-test)
  add_dependencies(docker-make-install ${DOCKER_IMAGE_AS_DIR}-make-install)
  add_dependencies(docker-make-clean ${DOCKER_IMAGE_AS_DIR}-make-clean)
  add_dependencies(docker-ctest ${DOCKER_IMAGE_AS_DIR}-ctest)

endmacro()
