#[=======================================================================[.rst:
ctest_tools.cmake
----
This file contains the functions used to setup ctest config

#]=======================================================================]

#[=======================================================================[.rst:
Call this function after each ctest step (ctest_configure,ctest_build ...)
to handle errors and submission to cdash

Usage :

post_ctest(PHASE <phase_name>)

with phase_name in (Configure, Build, Test).
#]=======================================================================]
function(post_ctest)
  set(oneValueArgs PHASE)
  set(options FORCE)
  cmake_parse_arguments(run "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  message("------> status/result : ${_STATUS}/${_RESULT}")
  # This function is supposed to submit to cdash
  # only if cmake terminates with an error (e.g. FATAL_ERROR) or ctest fails
  if(NOT _STATUS EQUAL 0 OR NOT _RESULT EQUAL 0 OR run_FORCE)
    if(CDASH_SUBMIT)
      ctest_submit(
        RETURN_VALUE RETURN_STATUS
        CAPTURE_CMAKE_ERROR SUBMISSION_STATUS
        )
      if(NOT SUBMISSION_STATUS EQUAL 0)
        message(WARNING " *** submission failure *** ")
      endif()
    else()
      message("- Results won't be submitted to a cdash server.\n")
      return()
    endif()
  endif()
  
  if(NOT _STATUS EQUAL 0 OR NOT _RESULT EQUAL 0)
    message(FATAL_ERROR "\n\n *** ${run_PHASE} process failed *** \n\n")
  endif()
  unset(_RESULT PARENT_SCOPE)
  unset(_STATUS PARENT_SCOPE)
  message("=============== End of ctest ${run_PHASE} =============== ")

endfunction()

# Set site name (cdash) according to current host status. 
function(set_site_name)
  
  # -- Query host system information --
  # --> to set ctest site for cdash.
  #include(cmake_host_system_information)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  cmake_host_system_information(RESULT fqdn QUERY FQDN)

  cmake_host_system_information(RESULT osname QUERY OS_NAME)
  cmake_host_system_information(RESULT osplatform QUERY OS_PLATFORM)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)

  string(STRIP ${osname} osname)
  string(STRIP ${osplatform} osplatform)

  if(CI_GITLAB)
    string(SUBSTRING $ENV{CI_JOB_IMAGE} 39 -1 dockerimagename) 
    string(STRIP ${dockerimagename} dockerimagename)
    set(hostname "[ ${dockerimagename} from gitlab registry]")
  endif()

  set(_SITE "${osname}-${osplatform}-host=${hostname}")

  string(STRIP _SITE ${_SITE})
  set(CTEST_SITE "${_SITE}" PARENT_SCOPE)
endfunction()

# set build name, according to host, ci, git status ...
function(set_cdash_build_name)
  # Get hash for commit of current version of Siconos
  # Saved by CI in CI_COMMIT_SHORT_SHA.
  if(CI_GITLAB)
    set(branch_commit "$ENV{CI_COMMIT_REF_NAME}-$ENV{CI_COMMIT_SHORT_SHA}")
  else()
    find_package(Git)
    execute_process(COMMAND
      ${GIT_EXECUTABLE} rev-parse --short HEAD
      OUTPUT_VARIABLE COMMIT_SHORT_SHA
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
    execute_process(COMMAND
      ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      OUTPUT_VARIABLE COMMIT_REF_NAME
      OUTPUT_STRIP_TRAILING_WHITESPACE
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
    set(branch_commit "${COMMIT_REF_NAME}-${COMMIT_SHORT_SHA}")
  endif()
  include(${CTEST_SOURCE_DIRECTORY}/cmake/SiconosVersion.cmake)  
  set(_name "Siconos(${SICONOS_VERSION}-devel,${branch_commit})")
  if(USER_OPTIONS_FILE)
    get_filename_component(_fname ${USER_OPTIONS_FILE} NAME)
    string(STRIP ${_fname} _fname)
    set(_name "${_name}-options=${_fname}")
    string(STRIP ${_name} _name)
  else() # Config = config_samples/default.cmake
    set(_name "${_name}-options=default.cmake")
    string(STRIP ${_name} _name)    
  endif()
  set(CTEST_BUILD_NAME "${_name}" PARENT_SCOPE)
endfunction()

# Write a note file for cdash server.
# Content :
# - info. regarding the runner, the system
# - siconos config (user option file)
function(write_notes)
  cmake_host_system_information(RESULT osrelease QUERY OS_RELEASE)
  cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
  file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt "- osrelease: ${osrelease}\n")

  if(CI_GITLAB)
     file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt "- Runner: ci-gitlab runner $ENV{CI_RUNNER_DESCRIPTION}\n")
  else()
    file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt "- host name: ${hostname}\n")
  endif() 
 
  file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt "\n\n ------- Siconos user options file ------\n\n")
  if(USER_OPTIONS_FILE)
    file(READ ${USER_OPTIONS_FILE} _options_file)
  else()
    file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt " Using default file (${CTEST_SOURCE_DIRECTORY}/config_samples/default.cmake):\n")
    file(READ ${CTEST_SOURCE_DIRECTORY}/config_samples/default.cmake _options_file)
  endif()
  file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt ${_options_file})
  file(APPEND ${CTEST_BINARY_DIRECTORY}/notes.txt " ------- End of Siconos user options file ------\n\n" )

  set(CTEST_NOTES_FILES ${CTEST_BINARY_DIRECTORY}/notes.txt PARENT_SCOPE)

endfunction()

