# ===========================================
#  Set some predefined compilation flags for
#  c compiler.
#
#  Usage :
#  include(CCompilerSetup)
# 
# ===========================================
include(cTools)
include(CheckSymbolExists)

detect_c_version(C_VERSION)


# ==== Warnings ===
# activate warnings
# and turn some of them to errors
if(DEV_MODE)
  # --- Clang ----
  #add_c_options("-Weverything" Clang) # like Wall and more
  #add_c_options("-Werror=unreachable-code" "Clang")
  add_c_options("-Wno-string-plus-int" "Clang")

  # --- All compilers but MSVC (Microsoft Visual C) ---
  if(NOT MSVC)
    add_c_options("-Wall")
  endif()
  add_c_options("-Wuninitialized")

  add_c_options("-Werror=unknown-warning-option" Clang)

  # --- Options for any compiler ----
  add_c_options("-Wextra -Wno-unused-parameter")
  # -- warnings to errors --
  add_c_options("-Werror=implicit-function-declaration")
  # should be supported only by Clang. The last statement is important, otherwise nothing compiles ...
  # MB: yes, nothing compiles
  if(DEV_MODE_STRICT)
    add_c_options("-Werror=conversion")
  endif()
  add_c_options("-Wno-sign-conversion")
  add_c_options("-Wno-error=sign-conversion")
  # ADD_C_OPTIONS("-Wno-error=shorten-64-to-32") # for clang
  add_c_options("-Werror=switch-bool")
  add_c_options("-Werror=logical-not-parentheses")
  add_c_options("-Werror=sizeof-array-argument")
  add_c_options("-Werror=bool-compare")
  add_c_options("-Werror=array-bounds")
  add_c_options("-Werror=format-invalid-specifier")
  add_c_options("-Werror=type-limits")
  add_c_options("-Werror=incompatible-pointer-types")
  add_c_options("-Werror=missing-prototypes")
  add_c_options("-Werror=return-type")
  add_c_options("-Wstrict-overflow=4")
  add_c_options("-Werror=strict-aliasing")
  add_c_options("-Werror=trampolines")
  add_c_options("-Werror=int-to-pointer-cast")
  add_c_options("-Werror=pointer-to-int-cast")

  if((NOT WITH_MECHANISMS) AND (NOT WITH_PYTHON_WRAPPER))
    add_c_options("-Werror=missing-declarations")
  endif()

  # ubuntu (at least) build with those
  add_c_options("-Wformat=2")
  add_c_options("-Werror=format-security")
endif()

# Compiler Specific
add_c_options("-diag-disable 654" "Intel")
if(NOT ICCOK)
  add_c_options("-D__aligned__=ignored" "Intel")
endif()

if(USE_SANITIZER MATCHES "asan")
  APPEND_C_FLAGS("-fsanitize=leak -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
elseif(USE_SANITIZER MATCHES "msan")
  APPEND_C_FLAGS("-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
elseif(USE_SANITIZER MATCHES "cfi")
  APPEND_C_FLAGS("-fsanitize=cfi -flto -fno-omit-frame-pointer -B ${CLANG_LD_HACK}")
endif()

# === Others options ===
if(C_VERSION STRLESS "201112L")
  set(C_STD_VERSION "c99")
else()
  # default C standart is c11 or newer
  set(C_STD_VERSION "c11")
endif()

if(NOT MSVC)
  add_c_options("-std=${C_STD_VERSION}")
  add_c_options("-x${C_STD_VERSION}")
endif()

if(LLVM_ANALYSE)
  set(CMAKE_C_OUTPUT_EXTENSION ".bc")
  set(CMAKE_C_FLAGS "-emit-llvm")
endif()

