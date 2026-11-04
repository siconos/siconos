# ================================================================
# Setup used to generate doc with CI
#
#
#
# ================================================================

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

# --- Build/compiling options ---
set(WARNINGS_LEVEL 0 CACHE STRING "Set compiler diagnostics level. 0: no warnings, 1: developer's minimal warnings, 2: strict level, warnings to errors and so on. Default =0")
option(WITH_TESTING "Enable 'make test' target" ON)

# --- Documentation setup ---
option(WITH_DOCUMENTATION "Build Documentation" ON)

# --- List of external libraries/dependencies to be searched (or not) ---
option(WITH_BULLET "compilation with Bullet Bindings. Default = OFF" ON)
option(WITH_OpenCASCADE "compilation with OpenCascade Bindings. Default = OFF" OFF)
option(WITH_MUMPS "Compilation with the MUMPS solver. Default = OFF" OFF)
option(WITH_UMFPACK "Compilation with the UMFPACK solver. Default = OFF" OFF)
option(WITH_SUPERLU "Compilation with the SuperLU solver. Default = OFF" OFF)
option(WITH_SUPERLU_MT "Compilation with the SuperLU solver, multithreaded version. Default = OFF" OFF)
option(WITH_FCLIB "link with fclib when this mode is enable. Default = ON" ON)
