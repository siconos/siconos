add_library(cxd SHARED cxd.cpp)

set_target_properties(cxd PROPERTIES SOVERSION 1)
target_compile_options(cxd PRIVATE -Wall)
target_compile_features(cxd PUBLIC cxx_std_20)

option(CXD_BUILD_TESTS "Build tests" ON)
if(CXD_BUILD_TESTS)
  add_subdirectory(test)
endif()

# The CUDA device backend (phase 5): opt-in, never a hard build dependency.
# The parity test compiles the device kernel header against the same
# committed golden corpus the core suite pins. nvcc 12.0 rejects gcc 13+ as
# its host compiler, so the host compiler is selectable; the spike validated
# clang-14 with this toolkit.
option(CXD_BUILD_DEVICE "Build the CUDA device backend parity test" OFF)
if(CXD_BUILD_DEVICE)
  set(CXD_CUDA_HOST_COMPILER "clang++-14"
      CACHE STRING "host compiler nvcc drives (the toolkit rejects gcc 13+)")
  set(CMAKE_CUDA_HOST_COMPILER "${CXD_CUDA_HOST_COMPILER}")
  enable_language(CUDA)
  add_executable(cxd_device_test device/cxd_device_test.cu)
  target_include_directories(cxd_device_test
      PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/src/test)
  target_compile_features(cxd_device_test PRIVATE cxx_std_17)
  target_compile_options(cxd_device_test PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
      --expt-relaxed-constexpr
      -fmad=false>)
  set_target_properties(cxd_device_test PROPERTIES
      CUDA_ARCHITECTURES 86 CUDA_STANDARD 17)
  target_compile_definitions(cxd_device_test
      PRIVATE CXD_TEST_GOLDEN_DIR="${CMAKE_SOURCE_DIR}/src/test/golden")
  add_test(t_device_golden_parity cxd_device_test)
endif()
