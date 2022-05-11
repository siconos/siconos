# - Fclib - 
if(WITH_FCLIB)
  # Two ways:
  # - User asks explicitely for a specific (already installed) version of fclib
  #   by providing fclib_ROOT on cmake command line
  # - else: use fetchcontent to download and install fclib as a siconos part
  message(STATUS "Start setup for fclib ... ")
  if(FCLIB_ROOT)
    find_package(FCLIB 3.0.0 CONFIG REQUIRED)
  else()
    message(STATUS "FCLIB will be downloaded from github repository and installed as a siconos component")
    FetchContent_Declare(fclib
      GIT_REPOSITORY    https://github.com/FrictionalContactLibrary/fclib.git
      GIT_TAG           origin/master
      GIT_SHALLOW TRUE
      UPDATE_DISCONNECTED TRUE # Do not update git repo at each run
      CMAKE_ARGS WITH_CXX ${WITH_CXX}
      LOG_CONFIGURE TRUE
      LOG_BUILD TRUE
      #     LOG_INSTALL TRUE
      )
    set(WITH_TESTS ${WITH_TESTING})
    FetchContent_MakeAvailable(fclib)
    add_library(FCLIB::fclib ALIAS fclib)
  endif()
  message(STATUS "End setup for fclib ... ")

  if(WITH_PYTHON_WRAPPER)
    include(swig_python_tools)
    set(COMPONENT fclib)
    #add_siconos_swig_sub_module("./fclib")
    set_property(SOURCE fclib.i PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)
  endif()

endif()

