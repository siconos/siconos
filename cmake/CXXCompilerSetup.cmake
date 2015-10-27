# ===========================================
#  Set some predefined compilation flags for
#  cxx compiler.
#
#  Usage :
#  include(CXXCompilerSetup)
# 
# ===========================================

if(NOT CMAKE_CXX_COMPILER OR NOT CMAKE_CXX_COMPILER_WORKS)
  message(ABORT "no cxx compiler")
endif()
  
include(cxxTools)

detect_cxx_standard(CXXVERSION)

# Check compiler version.
if(CMAKE_COMPILER_IS_GNUCXX)
  detect_cxx_compiler_version(compiler_version)
  if(compiler_version VERSION_LESS "4.4")
    message(FATAL_ERROR "gcc greater than 4.4 needed")
  endif()  
endif()

# ==== Warnings ===
# activate warnings
# and turn some of them to errors
if(DEV_MODE)
  # --- Clang ----
  #add_cxx_options("-Weverything" Clang) # like Wall and more
  add_cxx_options("-Werror=unreachable-code" "Clang")
  # --- All compilers but MSVC (Microsoft Visual C) ---
  if(NOT MSVC)
    add_cxx_options("-Wall")
  endif()
  add_cxx_options("-Wuninitialized")
  add_cxx_options("-Werror=unknown-warning-option" Clang)

  # --- Options for any compiler ----
  add_cxx_options("-Wextra -Wno-unused-parameter")
  # -- warnings to errors --
  add_cxx_options("-Werror=implicit-function-declaration")
  # should be supported only by Clang. The last statement is important, otherwise nothing compiles ...
  add_cxx_options("-Werror=conversion -Wno-sign-conversion -Wno-error=sign-conversion Wno-shorten-64-to-32 -Wno-error=shorten-64-to-32")
  # ADD_C_OPTIONS("-Wno-error=shorten-64-to-32") # for clang
  add_cxx_options("-Werror=switch-bool")
  add_cxx_options("-Werror=logical-not-parentheses")
  add_cxx_options("-Werror=sizeof-array-argument")
  add_cxx_options("-Werror=bool-compare")
  add_cxx_options("-Werror=array-bounds")
  add_cxx_options("-Werror=format-invalid-specifier")
  add_cxx_options("-Werror=type-limits")
  add_cxx_options("-Wodr")

  if((NOT WITH_MECHANISMS) AND (NOT WITH_PYTHON_WRAPPER))
    add_cxx_options("-Werror=missing-declarations")
  endif()
  if(NOT WITH_OCC AND NOT WITH_MECHANISMS)
    add_cxx_options("-Werror=overloaded-virtual")
  endif()
  add_cxx_options("-Wc++11-compat-deprecated-writable-strings")
endif()

# add_cxx_options("-static -static-libgcc -static-libstdc++" "GNU;Clang")
# way too verbose with MSVC

if(USE_SANITIZER MATCHES "asan")
  APPEND_CXX_FLAGS("-fsanitize=leak -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
elseif(USE_SANITIZER MATCHES "msan")
  APPEND_CXX_FLAGS("-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
endif()


# === Others options ===
if(NOT CXXVERSION STRLESS "201102L" AND DEV_MODE)
  add_cxx_options("-Wsuggest-final-types")
  add_cxx_options("-Wsuggest-final-methods")
  if(NOT WITH_PYTHON_WRAPPER)
    add_cxx_options("-Wzero-as-null-pointer-constant")
  endif()
endif()

# === Compiler Specific ===
# Intel
add_cxx_options("-diag-disable 654" "Intel")
if(NOT ICCOK)
  add_cxx_options("-D__aligned__=ignored" "Intel")
endif(NOT ICCOK)

# Clang
add_cxx_options("-Wno-string-plus-int" "Clang")

# Do we still need this?
append_cxx_flags("-D_NUMERICS_INTERNAL_CXX_")
