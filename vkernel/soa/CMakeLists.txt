#cmake_minimum_required(VERSION 3.2)
#project(
#  ${SKBUILD_PROJECT_NAME}
#  VERSION ${SKBUILD_PROJECT_VERSION}
#  LANGUAGES CXX C)

set(SOA_CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/vkernel/soa)

set(CMAKE_CUDA_FLAGS "-ccbin=${CMAKE_CXX_COMPILER} -allow-unsupported-compiler")

#set(WITH_GPU 1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckCXXCompilerFlag)
set(Siconos_DIR ${SICONOS_INSTALL_DIR})

find_package(Siconos REQUIRED)
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

include_directories(${SOA_CMAKE_SOURCE_DIR}/src)

if(WITH_GPU)
  include_directories(${SOA_CMAKE_SOURCE_DIR}/src/cuNSearch)
else()
  include_directories(${SOA_CMAKE_SOURCE_DIR}/src/CompactNSearch)
endif()

link_directories(${CMAKE_BUILD_DIR})

if(WITH_GPU)
  set(flags "-DWITH_GPU -DCUNSEARCH_USE_DOUBLE_PRECISION")
endif()

set(flags "${flags} -Wall --std=c++20")
set(fconcepts_diagnostics)
check_cxx_compiler_flag("-fbracket-depth=300" fbracket_depth)
check_cxx_compiler_flag("-fconcepts-diagnostics-depth=1" fconcepts_diagnostics)
check_cxx_compiler_flag("-ftemplate-backtrace-limit=1" ftemplate_backtrace_limit)
if(fconcepts_diagnostics)
  set(flags "${flags} -fconcepts-diagnostics-depth=1")
endif()
if(ftemplate_backtrace_limit)
  set(flags "${flags} -ftemplate-backtrace-limit=1")
endif()
if(fbracket_depth)
  set(flags "${flags} -fbracket-depth=300")
endif()

# check_cxx_compiler_flag("-fprofile-instr-generate" fprofile_instr_generate)
# check_cxx_compiler_flag("-fcoverage-mapping" fcoverage_mapping)

# if(fprofile_instr_generate)
#   set(flags "${flags} -fprofile-instr-generate -g")
# endif()
# if(fcoverage_mapping)
#   set(flags "${flags} -fcoverage-mapping")
# endif()

# python_add_library(_nothing MODULE ./src/siconos/pybind/nothing.cpp WITH_SOABI)
# set_target_properties(_nothing PROPERTIES COMPILE_FLAGS "${flags}")
# target_link_libraries(_nothing PRIVATE Eigen3::Eigen Siconos::numerics benchmark pthread fmt pybind11::headers)
# target_compile_definitions(_nothing PRIVATE VERSION_INFO=${PROJECT_VERSION})

# install(TARGETS _nothing DESTINATION nothing)

python_add_library(_nonos MODULE
  ./src/siconos/py/main.cpp
#  ./src/siconos/py/py_null.cpp
#  ./src/siconos/py/main_disks.cpp
#  ./src/siconos/py/py_attached_storages.cpp
#  ./src/siconos/py/py_attributes.cpp
#  ./src/siconos/py/py_add.cpp
#  ./src/siconos/py/py_methods.cpp
  WITH_SOABI)

set_target_properties(_nonos PROPERTIES COMPILE_FLAGS "${flags}")

if(WITH_GPU)
  target_link_libraries(_nonos PRIVATE mechanics Eigen3::Eigen numerics benchmark pthread fmt cuNSearch pybind11::headers)
else()
  target_link_libraries(_nonos PRIVATE mechanics Eigen3::Eigen numerics benchmark pthread fmt CompactNSearch pybind11::headers)
endif()


target_compile_definitions(_nonos PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS _nonos COMPONENT python DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos/)
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/py/nonos/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/src/siconos/py/nonos/bridge.py
  COMPONENT python
  DESTINATION ${SICONOS_PYTHON_INSTALL_DIR}/nonos/)

add_subdirectory(src)

# # use, i.e. don't skip the full RPATH for the build tree
# set(CMAKE_SKIP_BUILD_RPATH FALSE)

# # when building, don't use the install RPATH already
# # (but later on when installing)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


OPTION(BUILD_AS_SHARED_LIBS "Build all the libraries as shared" ON)

if(WITH_GPU)
  add_subdirectory(extern/cuNSearch)
else()
  add_subdirectory(extern/CompactNSearch)
endif()

