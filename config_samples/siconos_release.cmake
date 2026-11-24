# Release config based on siconos_asan.cmake
# --- List of siconos components to build and install ---
# Each component can be enabled/disabled with a WITH_<component> option.
# If the COMPONENTS variable is not set, it is derived from these options.
option(WITH_externals "Build externals component" ON)
option(WITH_numerics "Build numerics component" ON)
option(WITH_kernel "Build kernel component" ON)
option(WITH_control "Build control component" ON)
option(WITH_mechanics "Build mechanics component" ON)
option(WITH_mechanisms "Build mechanisms component" OFF)
option(WITH_io "Build io component" ON)
option(WITH_vkernel "Build vkernel component" OFF)
option(WITH_PYB11_WRAPPER "Build and install python bindings" ON)

set(WARNINGS_LEVEL 1 CACHE INTERNAL "Set compiler diagnostics level.")
option(WITH_TESTING "Enable 'make test' target" ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "build type (compiler options)" FORCE)

option(WITH_DOCUMENTATION "Build Documentation" OFF)

option(WITH_BULLET "compilation with Bullet Bindings" ON)
option(WITH_OpenCASCADE "compilation with OpenCascade Bindings" OFF)
option(WITH_MUMPS "Compilation with the MUMPS solver" OFF)
option(WITH_UMFPACK "Compilation with the UMFPACK solver" OFF)
option(WITH_SUPERLU "Compilation with the SuperLU solver" OFF)
option(WITH_SUPERLU_MT "Compilation with the SuperLU solver, multithreaded version" OFF)
option(WITH_FCLIB "link with fclib when this mode is enable" OFF)
option(WITH_SANITIZER "Activate asan" OFF)
