# Compilations and tests inside docker containers

# DOCKER_IMAGE specifications in ${CMAKE_SOURCE_DIR}/Build/Docker/${DOCKER_IMAGE}
# DOCKER_SHARED_DIRECTORIES : directories visibles in the container (source directory and binary directory are visibles)
# DOCKER_PROJECT_SOURCE_DIR : project source directory (for example some git clone)
# DOCKER_CMAKE_FLAGS : flags passed to cmake 
# DOCKER_MAKE_FLAGS  : flags passed to make
# DOCKER_MAKE_TEST_FLAGS  : flags passed to make test

macro(add_docker_targets)

  set(options)
  set(oneValueArgs DOCKER_TEMPLATES DOCKER_TEMPLATE DOCKER_IMAGE DOCKER_IMAGE_DIR DOCKER_REPOSITORY DOCKER_FILE DOCKER_WORKDIR DOCKER_HOST_INSTALL_PREFIX DOCKER_SHARED_DIRECTORIES DOCKER_CMAKE_FLAGS DOCKER_MAKE_FLAGS DOCKER_MAKE_TEST_FLAGS DOCKER_CTEST_DRIVER DOCKER_HOSTNAME DOCKER_COMMAND DOCKER_CMAKE_WRAPPER DOCKER_PROJECT_SOURCE_DIR)

  string(REPLACE ":" "-" DOCKER_DISTRIB_AS_NAME ${DOCKER_DISTRIB})
  string(REPLACE "." "-" DOCKER_DISTRIB_AS_NAME ${DOCKER_DISTRIB_AS_NAME})

  if(DOCKER_TEMPLATES)
    string(REPLACE "," ";" DOCKER_TEMPLATES_LIST ${DOCKER_TEMPLATES})

    string(REPLACE "," "-" DOCKER_IMAGE ${DOCKER_TEMPLATES})
    string(REPLACE "+" "x" DOCKER_IMAGE ${DOCKER_IMAGE})

    set(DOCKER_IMAGE ${DOCKER_DISTRIB_AS_NAME}-${DOCKER_IMAGE})
  endif()    

  if(DOCKER_TEMPLATE)
    set(DOCKER_IMAGE ${DOCKER_DISTRIB_AS_NAME}-${DOCKER_TEMPLATE})
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
    set(DOCKER_IMAGE_DIR .)
  endif()
  
  if(NOT DOCKER_REPOSITORY)
    message(FATAL_ERROR "Docker : DOCKER_REPOSITORY unset")
  endif()

  # Default behavior: use all components 
  if(NOT SICONOS_COMPONENTS)
	include(Tools)
	set_components(externals;numerics;kernel;control;mechanics;io)
  endif()
  # ================== Create Dockerfile ==================
  # 
  # Set Dockerfile name
  set(GENERATED_DOCKER_FILE ${CMAKE_CURRENT_BINARY_DIR}/Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR}/Dockerfile)

  file(REMOVE ${GENERATED_DOCKER_FILE})
  
  message(STATUS "Docker templates list is : ${DOCKER_TEMPLATES_LIST}")

  # Check for templates
  foreach(_dt ${DOCKER_TEMPLATES_LIST})
    if(EXISTS ${DOCKER_IMAGE_DIR}/${_dt})
      file(READ ${DOCKER_IMAGE_DIR}/${_dt} _contents_)
      file(APPEND ${GENERATED_DOCKER_FILE} ${_contents_})
      list(REMOVE_ITEM DOCKER_TEMPLATES_LIST ${_dt})
    endif()
  endforeach()

  if(NOT DOCKER_MKSENV_SPLIT)
    set(DOCKER_MKSENV_SPLIT false)
  endif()

  # Call print_command from mksenv.py to write Dockerfile
  execute_process(COMMAND ${DOCKER_MKSENV_COMMAND} --docker --distrib ${DOCKER_DISTRIB} --pkgs ${DOCKER_TEMPLATES} --split=${DOCKER_MKSENV_SPLIT} ${DOCKER_MKSENV_INPUT} OUTPUT_VARIABLE _contents_)
  file(APPEND ${GENERATED_DOCKER_FILE} ${_contents_})
  file(APPEND ${GENERATED_DOCKER_FILE} "RUN mkdir -p /usr/local\n")
  # ================== End of Dockerfile creation ==================
  

  if(NOT DOCKER_WORKDIR)
    set(DOCKER_WORKDIR /docker-workdir)
  endif()
  
  set(DOCKER_WORKDIR_VOLUME "workdir-${DOCKER_PROJECT_SOURCE_DIR}-${DOCKER_IMAGE_AS_DIR}")
  string(REPLACE "/" "-" DOCKER_WORKDIR_VOLUME ${DOCKER_WORKDIR_VOLUME})
  string(REPLACE ":" "-" DOCKER_WORKDIR_VOLUME ${DOCKER_WORKDIR_VOLUME})
  string(REPLACE "." "-" DOCKER_WORKDIR_VOLUME ${DOCKER_WORKDIR_VOLUME})

  set(CTEST_OPTIONS_WITHOUT_DOCKER)
  # List of options (-DOPTION=VALUE ...) to be passed to ctest
  # Warning : these are not siconos conf. options for cmake but
  # the options needed by ctest to prepare the pipeline (cmake, make ...)
  # like SITE, config file ...
  # If CI_OPTIONS is empty, default.cmake will be used.
  set(CTEST_OPTIONS ${CI_OPTIONS};-DWITH_DOCKER=0;-V)


  foreach(_f ${CTEST_OPTIONS})
	string(REGEX MATCH "^-DDOCKER_.*" _mf ${_f})
	if (NOT _mf)
	  list(APPEND CTEST_OPTIONS_WITHOUT_DOCKER ${_f})
	endif()
  endforeach()

  
  set(DOCKER_VFLAGS)
  foreach(_D ${DOCKER_SHARED_DIRECTORIES};${DOCKER_PROJECT_SOURCE_DIR})
    string(REGEX MATCH ":" _md ${_D})
    if(_md)
      set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${_D})
    else()
      set(DOCKER_VFLAGS ${DOCKER_VFLAGS} -v ${_D}:${_D})
    endif()
  endforeach()

  if(NOT DOCKER_HOSTNAME)
    set(DOCKER_HOSTNAME ${DOCKER_IMAGE_AS_DIR})
    string(LENGTH ${DOCKER_HOSTNAME} DOCKER_HOSTNAME_LENGTH)
    if(${DOCKER_HOSTNAME_LENGTH} GREATER 63)
      string(SUBSTRING ${DOCKER_HOSTNAME} 0 63 DOCKER_HOSTNAME)
    endif()
  endif()

  IF(DOCKER_CMAKE_WRAPPER)
    message(STATUS "Docker cmake command : ${DOCKER_CMAKE_WRAPPER}")
  ELSE(DOCKER_CMAKE_WRAPPER)
    SET(DOCKER_CMAKE_WRAPPER "cmake")
  ENDIF(DOCKER_CMAKE_WRAPPER)

  message(STATUS "ctest flags : ${CTEST_OPTIONS_WITHOUT_DOCKER}")
  message(STATUS "Docker make flags : ${DOCKER_MAKE_FLAGS}")
  message(STATUS "Docker make test flags : ${DOCKER_MAKE_TEST_FLAGS}")
  message(STATUS "Docker hostname : ${DOCKER_HOSTNAME}")


  # ================== Create targets ==================

  # 
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-clean-usr-local
    COMMENT "Docker clean usr-local : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} rm -f -v ${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local || /bin/true
    )

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-clean
    COMMENT "Docker clean : ${DOCKER_IMAGE}"
    DEPENDS ${DOCKER_IMAGE_AS_DIR}-clean-usr-local
    COMMAND ${DOCKER_COMMAND} rm -f -v ${DOCKER_WORKDIR_VOLUME} || /bin/true
    COMMAND ${DOCKER_COMMAND} rmi -f ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} || /bin/true
    )

  # == Target to build a docker image and create a container ==
  # - working dir = Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR}
  # - Dockerfile path =  ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} (i.e. use GENERATED_DOCKER_FILE created above)
  # make docker-build will:
  #  * build a new docker image
  #  * create /usr/local and DOCKER_WORKDIR volumes.
  # Those volumes will be mounted when running image/container and will be persistent
  # image ...
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-build
    COMMENT "Docker Build : ${DOCKER_IMAGE}"
    COMMAND cd Docker/Context/${DOCKER_REPOSITORY}/${DOCKER_IMAGE_AS_DIR} && docker build -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} .
    )

  # bind DOCKER_WORKDIR inside container
  add_custom_command(
    TARGET ${DOCKER_IMAGE_AS_DIR}-build
    PRE_BUILD 
    COMMENT "docker create workdir"
    COMMAND ${DOCKER_COMMAND} create --name=${DOCKER_WORKDIR_VOLUME} -v ${DOCKER_WORKDIR} -i -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} true && echo done || echo already done)

  # bind usr/local inside container
  add_custom_command(
    TARGET ${DOCKER_IMAGE_AS_DIR}-build
    PRE_BUILD 
    COMMENT "docker create /usr/local"
    COMMAND ${DOCKER_COMMAND} create --name=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local -v /usr/local -i -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} true && echo done || echo already done
    )

  # == target to run cmake inside previously created docker container ==
  # Configure siconos inside docker container.
  # Note : --rm=true is supposed to cleanup container file system after exit. But
  # volumes /usr/local and WORKDIR are persistent
  # (see docker doc, https://docs.docker.com/engine/reference/run/#clean-up---rm,
  # : "if the original volume was specified with a name it will not be removed."
  # 
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-cmake
    COMMENT "Docker cmake : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} ${DOCKER_CMAKE_WRAPPER} ${DOCKER_PROJECT_SOURCE_DIR} ${SICONOS_CMAKE_OPTIONS} -DCOMPONENTS=${SICONOS_COMPONENTS} VERBATIM)
  
  # == target to run make inside previously created docker container ==
  # Build siconos inside docker container
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make
    COMMENT "Docker make : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_FLAGS})
  
  # == target to run make test inside previously created docker container ==
  # Make test (siconos) inside docker container
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-test
    COMMENT "Docker make test : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_TEST_FLAGS} test)
  
  # == target to install siconos inside previously created docker container ==
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-install
    COMMENT "Docker make install : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_INSTALL_FLAGS} -ki install)

  # == target to uninstall siconos inside previously created docker container ==
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-uninstall
    COMMENT "Docker make uninstall : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_INSTALL_FLAGS} -ki uninstall)

  # == target to create doc inside previously created docker container ==
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-doc
    COMMENT "Docker make doc : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_DOC_FLAGS} -ki doc)

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-upload
    COMMENT "Docker make upload : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_UPLOAD_FLAGS} -ki upload)
  
  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-make-clean
    COMMENT "Docker make clean : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} make ${DOCKER_MAKE_CLEAN_FLAGS} clean)

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-ctest
    COMMENT "Docker ctest : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} ctest -DCTEST_SOURCE_DIRECTORY=${DOCKER_PROJECT_SOURCE_DIR} -DCTEST_BINARY_DIRECTORY=${DOCKER_WORKDIR} -S ${DOCKER_CTEST_DRIVER} -DSITE=${DOCKER_HOSTNAME} -DCMAKE_WRAPPER=${DOCKER_CMAKE_WRAPPER} ${CTEST_OPTIONS_WITHOUT_DOCKER})

  add_custom_target(
    ${DOCKER_IMAGE_AS_DIR}-interactive
    COMMENT "Docker interactive : ${DOCKER_IMAGE}"
    COMMAND ${DOCKER_COMMAND} run -h ${DOCKER_HOSTNAME} --rm=true ${DOCKER_VFLAGS} --volumes-from=${DOCKER_WORKDIR_VOLUME} --volumes-from=${DOCKER_REPOSITORY}-${DOCKER_IMAGE}-usr-local --workdir=${DOCKER_WORKDIR} -i -t ${DOCKER_REPOSITORY}/${DOCKER_IMAGE} /bin/bash)

  add_custom_target(
    docker-hard-clean
    COMMENT "Docker cleanup (containers, unused images and volumes)"
    COMMAND sh ${CMAKE_SOURCE_DIR}/docker-cleanup.sh
    )

  if(NOT TARGET docker-clean-usr-local)
    add_custom_target(
      docker-clean-usr-local
      COMMENT "Docker clean"
      )
  endif()

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

  if(NOT TARGET docker-make-doc)
    add_custom_target(
      docker-make-doc
      COMMENT "Docker make doc"
      )
  endif()

  if(NOT TARGET docker-make-upload)
    add_custom_target(
      docker-make-upload
      COMMENT "Docker make upload"
      )
  endif()

  if(NOT TARGET docker-make-clean)
    add_custom_target(
      docker-make-clean
      COMMENT "Docker make clean"
      )
  endif()

  if(NOT TARGET docker-make-uninstall)
	add_custom_target(
	  docker-make-uninstall
	  COMMENT "Docker make uninstall"
	  )
  endif()
  
  if(NOT TARGET docker-ctest)
    add_custom_target(
      docker-ctest
      COMMENT "Docker ctest"
      )
  endif()

  if(NOT TARGET docker-interactive)
    add_custom_target(
      docker-interactive
      COMMENT "Docker interactive"
      )
  endif()

  add_dependencies(docker-clean-usr-local ${DOCKER_IMAGE_AS_DIR}-clean-usr-local)
  add_dependencies(docker-clean ${DOCKER_IMAGE_AS_DIR}-clean)
  add_dependencies(docker-build ${DOCKER_IMAGE_AS_DIR}-build)
  add_dependencies(docker-cmake ${DOCKER_IMAGE_AS_DIR}-cmake)
  add_dependencies(docker-make ${DOCKER_IMAGE_AS_DIR}-make)
  add_dependencies(docker-make-test ${DOCKER_IMAGE_AS_DIR}-make-test)
  add_dependencies(docker-make-install ${DOCKER_IMAGE_AS_DIR}-make-install)
  add_dependencies(docker-make-uninstall ${DOCKER_IMAGE_AS_DIR}-make-uninstall)
  add_dependencies(docker-make-doc ${DOCKER_IMAGE_AS_DIR}-make-doc)
  add_dependencies(docker-make-upload ${DOCKER_IMAGE_AS_DIR}-make-upload)
  add_dependencies(docker-make-clean ${DOCKER_IMAGE_AS_DIR}-make-clean)
  add_dependencies(docker-ctest ${DOCKER_IMAGE_AS_DIR}-ctest)
  add_dependencies(docker-interactive ${DOCKER_IMAGE_AS_DIR}-interactive)

endmacro()
