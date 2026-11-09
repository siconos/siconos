include(CheckCXXCompilerFlag)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
link_directories(${CMAKE_BUILD_DIR})

add_executable(tground tground.cpp)
add_executable(tpattern tpattern.cpp)

set(flags "-Wall --std=c++20")
set(fconcepts_diagnostics)
check_cxx_compiler_flag("-fconcepts-diagnostics-depth=0" fconcepts_diagnostics)
if(fconcepts_diagnostics)
  set(flags "${flags} -fconcepts-diagnostics-depth=0")
endif()

set_target_properties(tground PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
set_target_properties(tpattern PROPERTIES COMPILE_FLAGS "${flags}")# -fconcepts-diagnostics-depth=10")
