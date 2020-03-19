# ===========================================
#  Set some predefined compilation flags for
#  cxx compiler.
#
#  Usage :
#
#  project(... CXX ...)
#  include(CXXCompilerSetup)
# 
# ===========================================

include(cxxTools)

# Version of C++
# detect_cxx_standard(CXXVERSION)

# Used to avoid tests below if same as previous run
#set("CXXVERSION_CUR" "${CXXVERSION}:${CMAKE_CXX_STANDARD}")

# Used to avoid tests above if same as previous run
#set("CXXVERSION_LAST" "${CXXVERSION_CUR}"
#    CACHE INTERNAL "Last value of CXXVERSION" FORCE)

# Check compiler version.
if(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.5")
    message(FATAL_ERROR "gcc greater than 4.5 needed")
  endif()  
endif()

# ==== Warnings ===
# activate warnings
# and turn some of them to errors
if(DEV_MODE)
  # --- Clang ----
  #add_cxx_options("-Weverything" Clang) # like Wall and more
  if (WITH_BULLET)
    # some unreachable code errors come from the Bullet headers
    add_cxx_options("-Wunreachable-code" "Clang")
  else()
    add_cxx_options("-Werror=unreachable-code" "Clang")
  endif()
  # --- All compilers but MSVC (Microsoft Visual C) ---
  if(NOT MSVC)
    add_cxx_options("-Wall")
  endif()
  # add_cxx_options("-Wuninitialized") # Done by wall
  add_cxx_options("-Werror=unknown-warning-option" Clang)
  add_cxx_options("-Wno-unused-local-typedef")  # remove ? 
  # --- Options for any compiler ----
  add_cxx_options("-Wextra -Wno-unused-parameter")
  # -- warnings to errors --

  # should be supported only by Clang. The last statement is important, otherwise nothing compiles ...
  # MB: yes, nothing compiles
  if(DEV_MODE_STRICT)
    add_cxx_options("-Werror=conversion")
  endif()
  add_cxx_options("-Wno-sign-conversion")
  add_cxx_options("-Wno-error=sign-conversion")
  add_cxx_options("-Wno-shorten-64-to-32") # unknown by gcc
  add_cxx_options("-Wno-error=shorten-64-to-32") # unknown by gcc
  add_cxx_options("-Werror=switch-bool")
  add_cxx_options("-Werror=logical-not-parentheses")
  add_cxx_options("-Werror=sizeof-array-argument")
  add_cxx_options("-Werror=bool-compare")
  add_cxx_options("-Werror=array-bounds")
  add_cxx_options("-Werror=format-invalid-specifier") # unknown by gcc
  add_cxx_options("-Werror=type-limits")
  add_cxx_options("-Werror=return-type")
  # add_cxx_options("-Wodr") # enable by default, useless
  if((NOT WITH_OCE) AND (NOT WITH_PYTHON_WRAPPER))
    add_cxx_options("-Werror=missing-declarations")
  endif()
  if(NOT WITH_OCE AND NOT WITH_SERIALIZATION)
    add_cxx_options("-Werror=overloaded-virtual")
  endif()

  add_cxx_options("-Wc++11-compat-deprecated-writable-strings") # obsolete ?
  # ubuntu (at least) build with those
  add_cxx_options("-Wformat=2")
  add_cxx_options("-Werror=format-security")

  if(NOT WITH_OCE)
    add_cxx_options("-Werror=non-virtual-dtor")
  endif()
endif()

IF(USE_LIBCXX) # Used for msan ?? 
  APPEND_CXX_FLAGS("-stdlib=libc++ -I${USE_LIBCXX}/include -I${USE_LIBCXX}/include/c++/v1")
  SET(_LIBCXX_FLAGS_TO_ADD "-L${USE_LIBCXX}/lib -lc++abi -Wl,-rpath,${USE_LIBCXX}/lib")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
ENDIF()


# === Others options ===
add_cxx_options("-Wsuggest-final-types")
add_cxx_options("-Wsuggest-final-methods")
if(NOT WITH_PYTHON_WRAPPER)
  add_cxx_options("-Wzero-as-null-pointer-constant")
endif()

# === Compiler Specific ===

# Clang
add_cxx_options("-Wno-string-plus-int" "Clang")

if(WITH_SERIALIZATION)
  add_cxx_options("-ftemplate-depth=1024" Clang)
endif(WITH_SERIALIZATION)
message(" CXX Flags (cmake setup) : ${CMAKE_CXX_FLAGS}")
