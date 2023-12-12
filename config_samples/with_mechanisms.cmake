# ===============================================================================
# Siconos with mechanisms component (requires OpenCascade !!)
#
#
# ===============================================================================

# --- List of siconos components to build and install ---
set(COMPONENTS externals numerics kernel control mechanics mechanisms io CACHE INTERNAL "List of siconos components to build and install")
option(WITH_PYTHON_WRAPPER "Build and install python bindings" ON)

# --- Build/compiling options ---
set(WARNINGS_LEVEL 0 CACHE INTERNAL "Set compiler diagnostics level. 0: no warnings, 1: developer's minimal warnings, 2: strict level, warnings to errors and so on.")
option(WITH_TESTING "Enable 'make test' target" ON)

# --- Documentation setup ---
option(WITH_DOCUMENTATION "Build Documentation" OFF)

# --- List of external libraries/dependencies to be searched (or not) ---
option(WITH_BULLET "compilation with Bullet Bindings" ON)
option(WITH_OCC "compilation with OCE Bindings" OFF)
option(WITH_MUMPS "Compilation with the MUMPS solver" OFF)
option(WITH_UMFPACK "Compilation with the UMFPACK solver" OFF)
option(WITH_SUPERLU "Compilation with the SuperLU solver" OFF)
option(WITH_SUPERLU_MT "Compilation with the SuperLU solver, multithreaded version" OFF)
option(WITH_FCLIB "link with fclib when this mode is enable" OFF)


