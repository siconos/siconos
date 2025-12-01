# - Fclib - 
if(WITH_FCLIB)
  # Three ways:
  # - No info about fclib location, try standard paths or use fetchcontent
  # - User asks explicitely for a specific (already installed) version of fclib
  #   by providing fclib_ROOT on cmake command line
  # - else: use fetchcontent to download and install fclib as a siconos part
  find_package(fclib 3.1.0 CONFIG)

  if(NOT fclib_FOUND)
    message(STATUS "  fclib will be downloaded from github repository and installed as a siconos component")
    FetchContent_Declare(fclib
      GIT_REPOSITORY    https://github.com/FrictionalContactLibrary/fclib.git
      GIT_TAG           origin/master
      GIT_SHALLOW TRUE
      UPDATE_DISCONNECTED TRUE # Do not update git repo at each run
      CMAKE_ARGS 
        -DWITH_CXX=${WITH_CXX}
        -DFCLIB_WITH_MERIT_FUNCTIONS=${TRUC_WITH_MERIT_FUNCTIONS}
        -DWITH_TESTS=${WITH_TESTING}
      LOG_CONFIGURE TRUE
      LOG_BUILD TRUE
      #     LOG_INSTALL TRUE
      )
    set(WITH_TESTS ${WITH_TESTING})
    FetchContent_MakeAvailable(fclib)
    add_library(fclib::fclib ALIAS fclib)
    target_compile_options(fclib INTERFACE
    $<$<CXX_COMPILER_ID:GNU>:-w>
    $<$<CXX_COMPILER_ID:Clang>:-w>
    )
  endif()
  if(NOT TARGET fclib::fclib)
    add_library(fclib::fclib ALIAS fclib)
  endif()

  if(WITH_PYTHON_WRAPPER)
    include(swig_python_tools)
    set(COMPONENT fclib)
    #add_siconos_swig_sub_module("./fclib")
    set_property(SOURCE fclib.i PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)
  endif()

endif()

