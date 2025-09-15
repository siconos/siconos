# ================================================================
# Setup used to generate doc with CI
#
#
#
# ================================================================

# --- List of siconos components to build and install ---
set(COMPONENTS externals numerics kernel control mechanics io CACHE INTERNAL "List of siconos components to build and install")
option(WITH_PYTHON_WRAPPER "Build and install python bindings using swig. Default = ON" ON)

# --- Build/compiling options ---
set(WARNINGS_LEVEL 0 CACHE INTERNAL "Set compiler diagnostics level. 0: no warnings, 1: developer's minimal warnings, 2: strict level, warnings to errors and so on. Default =0")
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

