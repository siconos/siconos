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
#add_definitions(-DBOOST_NOEXCEPT)
# For SiconosConfig.h
option(SICONOS_USE_BOOST_FOR_CXX11 "Prefer BOOST features over C++ standard features even if C++xy is enabled" ON)
option(SICONOS_USE_MAP_FOR_HASH "Prefer std::map to std::unordered_map even if C++xy is enabled" ON)

include(cxxTools)

# Version of C++
detect_cxx_standard(CXXVERSION)

# Used to avoid tests below if same as previous run
set("CXXVERSION_CUR" "${CXXVERSION}:${CMAKE_CXX_STANDARD}")

# For SiconosConfig.h
detect_cxx_std_shared_ptr(SICONOS_STD_SHARED_PTR)

detect_cxx_std_array(SICONOS_STD_ARRAY)
detect_cxx_std_unordered_map(SICONOS_STD_UNORDERED_MAP)
detect_cxx_std_tuple(SICONOS_STD_TUPLE)
detect_cxx_std_to_string(SICONOS_STD_TO_STRING)
detect_cxx_std_functional(SICONOS_STD_FUNCTIONAL)

# Used to avoid tests above if same as previous run
set("CXXVERSION_LAST" "${CXXVERSION_CUR}"
    CACHE INTERNAL "Last value of CXXVERSION" FORCE)

# Check compiler version.
if(CMAKE_COMPILER_IS_GNUCXX)
  detect_cxx_compiler_version(compiler_version)
  if(compiler_version VERSION_LESS "4.5")
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
  add_cxx_options("-Wuninitialized")
  add_cxx_options("-Werror=unknown-warning-option" Clang)
  add_cxx_options("-Wno-unused-local-typedef")
  # --- Options for any compiler ----
  add_cxx_options("-Wextra -Wno-unused-parameter")
  # -- warnings to errors --
  add_cxx_options("-Werror=implicit-function-declaration")
  # should be supported only by Clang. The last statement is important, otherwise nothing compiles ...
  # MB: yes, nothing compiles
  if(DEV_MODE_STRICT)
    add_cxx_options("-Werror=conversion")
  endif()
  add_cxx_options("-Wno-sign-conversion")
  add_cxx_options("-Wno-error=sign-conversion")
  add_cxx_options("-Wno-shorten-64-to-32")
  add_cxx_options("-Wno-error=shorten-64-to-32")
  # ADD_C_OPTIONS("-Wno-error=shorten-64-to-32") # for clang
  add_cxx_options("-Werror=switch-bool")
  add_cxx_options("-Werror=logical-not-parentheses")
  add_cxx_options("-Werror=sizeof-array-argument")
  add_cxx_options("-Werror=bool-compare")
  add_cxx_options("-Werror=array-bounds")
  add_cxx_options("-Werror=format-invalid-specifier")
  add_cxx_options("-Werror=type-limits")
  add_cxx_options("-Werror=return-type")
  add_cxx_options("-Wodr")

  if((NOT WITH_MECHANISMS) AND (NOT WITH_PYTHON_WRAPPER))
    add_cxx_options("-Werror=missing-declarations")
  endif()
  if(NOT WITH_OCC AND NOT WITH_MECHANISMS AND NOT WITH_SERIALIZATION)
    add_cxx_options("-Werror=overloaded-virtual")
  endif()

  add_cxx_options("-Wc++11-compat-deprecated-writable-strings")
  # ubuntu (at least) build with those
  add_cxx_options("-Wformat=2")
  add_cxx_options("-Werror=format-security")

  if(NOT WITH_OCC AND NOT WITH_MECHANISMS)
    add_cxx_options("-Werror=non-virtual-dtor")
  endif()
endif()

# add_cxx_options("-static -static-libgcc -static-libstdc++" "GNU;Clang")
# way too verbose with MSVC

if(USE_SANITIZER MATCHES "asan")
  APPEND_CXX_FLAGS("-fsanitize=leak -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
elseif(USE_SANITIZER MATCHES "msan")
  APPEND_CXX_FLAGS("-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
elseif(USE_SANITIZER MATCHES "cfi")
  APPEND_CXX_FLAGS("-fsanitize=cfi -flto -fno-omit-frame-pointer -B ${CLANG_LD_HACK}")
endif()

IF(USE_LIBCXX)
  APPEND_CXX_FLAGS("-stdlib=libc++ -I${USE_LIBCXX}/include -I${USE_LIBCXX}/include/c++/v1")
  SET(_LIBCXX_FLAGS_TO_ADD "-L${USE_LIBCXX}/lib -lc++abi -Wl,-rpath,${USE_LIBCXX}/lib")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
ENDIF(USE_LIBCXX)


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

if(WITH_SERIALIZATION)
  add_cxx_options("-ftemplate-depth=1024" Clang)
endif(WITH_SERIALIZATION)
