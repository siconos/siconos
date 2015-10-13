# apply misc flags

IF(CMAKE_CXX_COMPILER_WORKS)
  INCLUDE(cxxVersion)
ENDIF()
INCLUDE(cVersion)

macro(ADD_CXX_OPTIONS OPT)

 STRING(REGEX REPLACE " " "" OPT_SANE "${OPT}")
 CHECK_CXX_ACCEPTS_FLAG("${OPT} ${_EXTRA_WARNING_FLAGS}" CXX_HAVE_${OPT_SANE})

 set(_compilers ${ARGN})
 IF(_compilers)
  SET(ADD_OPTION FALSE)
  FOREACH(_compiler ${_compilers})
   IF (${CMAKE_CXX_COMPILER_ID} MATCHES ${_compiler})
    SET(ADD_OPTION TRUE)
   ENDIF()
  ENDFOREACH()
 ELSE(_compilers)
  SET(ADD_OPTION TRUE)
 ENDIF(_compilers)

 IF(ADD_OPTION AND CXX_HAVE_${OPT_SANE})
  APPEND_CXX_FLAGS("${OPT}")
 ENDIF(ADD_OPTION AND CXX_HAVE_${OPT_SANE})

endmacro(ADD_CXX_OPTIONS)

macro(ADD_C_OPTIONS OPT)

 STRING(REGEX REPLACE " " "" OPT_SANE "${OPT}")
 CHECK_C_COMPILER_FLAG("${OPT} ${_EXTRA_WARNING_FLAGS}" C_HAVE_${OPT_SANE})

 set(_compilers ${ARGN})
 IF(_compilers)
  SET(ADD_OPTION FALSE)
  FOREACH(_compiler ${_compilers})
   IF (${CMAKE_C_COMPILER_ID} MATCHES ${_compiler})
    MESSAGE(STATUS "Adding option for compiler ${_compiler}")
    SET(ADD_OPTION TRUE)
   ENDIF()
  ENDFOREACH()
 ELSE(_compilers)
  SET(ADD_OPTION TRUE)
 ENDIF(_compilers)

 IF(ADD_OPTION AND C_HAVE_${OPT_SANE})
  APPEND_C_FLAGS("${OPT}")
 ENDIF(ADD_OPTION AND C_HAVE_${OPT_SANE})

endmacro(ADD_C_OPTIONS)


IF(CMAKE_C_COMPILER)
 INCLUDE(CheckCCompilerFlag)

 IF(${CMAKE_C_COMPILER_ID} MATCHES "Clang")
  SET(_EXTRA_WARNING_FLAGS "-Werror=unknown-warning-option")
 ELSE()
  SET(_EXTRA_WARNING_FLAGS "")
 ENDIF()

 detect_c_version(C_VERSION)

 IF(C_VERSION STRLESS "201112L")
  SET(C_STD_VERSION "c99")
 ELSE(C_VERSION STRLESS "201112L")
  # default C standart is c11 or newer
  SET(C_STD_VERSION "c11")
 ENDIF(C_VERSION STRLESS "201112L")

 IF(NOT MSVC)
  ADD_C_OPTIONS("-std=${C_STD_VERSION}")
  ADD_C_OPTIONS("-x${C_STD_VERSION}")
 ENDIF(NOT MSVC)

 # ADD_C_OPTIONS("-static -static-libgcc" "GNU;Clang")
 # way too verbose with MSVC
 IF(NOT MSVC)
  ADD_C_OPTIONS("-Wall")
 ENDIF(NOT MSVC)
 ADD_C_OPTIONS("-Werror=overloaded-virtual")
 ADD_C_OPTIONS("-Wextra -Wno-unused-parameter")
 ADD_C_OPTIONS("-Werror=implicit-function-declaration")
 ADD_C_OPTIONS("-Werror=conversion -Wno-sign-conversion -Wno-error=sign-conversion")
 # ADD_C_OPTIONS("-Wno-error=shorten-64-to-32") # for clang
 ADD_C_OPTIONS("-Werror=switch-bool")
 ADD_C_OPTIONS("-Werror=logical-not-parentheses")
 ADD_C_OPTIONS("-Werror=sizeof-array-argument")
 ADD_C_OPTIONS("-Werror=bool-compare")
 ADD_C_OPTIONS("-Werror=array-bounds")
 ADD_C_OPTIONS("-Werror=format-invalid-specifier")
 ADD_C_OPTIONS("-Werror=type-limits")
 ADD_C_OPTIONS("-Werror=incompatible-pointer-types")
 ADD_C_OPTIONS("-Werror=empty-body")
 ADD_C_OPTIONS("-Werror=implicit")
 ADD_C_OPTIONS("-Werror=pointer-to-int-cast")
 ADD_C_OPTIONS("-Werror=int-to-pointer-cast")
 # C specific
 ADD_C_OPTIONS("-Werror=missing-prototypes")

 # Compiler Specific
 ADD_C_OPTIONS("-diag-disable 654" "Intel")
 IF(NOT ICCOK)
  ADD_C_OPTIONS("-D__aligned__=ignored" "Intel")
 ENDIF(NOT ICCOK)

 ADD_C_OPTIONS("-Wno-string-plus-int" "Clang")
 ADD_C_OPTIONS("-Werror=unreachable-code" "Clang")

 IF(USE_SANITIZER MATCHES "asan")
   APPEND_C_FLAGS("-fsanitize=leak -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
 ELSEIF(USE_SANITIZER MATCHES "msan")
   APPEND_C_FLAGS("-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
 ENDIF(USE_SANITIZER MATCHES "asan")

ENDIF(CMAKE_C_COMPILER)


IF(CMAKE_CXX_COMPILER_WORKS)
 
  if (NOT CMAKE_CXX_COMPILER)
    message(ABORT "no cxx compiler")
  endif()
  
  detect_cxx_version(CXXVERSION)

  INCLUDE(TestCXXAcceptsFlag)
  
  IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    SET(_EXTRA_WARNING_FLAGS "-Werror=unknown-warning-option")
  ELSE()
    SET(_EXTRA_WARNING_FLAGS "")
  ENDIF()

  # ADD_CXX_OPTIONS("-static -static-libgcc -static-libstdc++" "GNU;Clang")
  # way too verbose with MSVC
  IF(NOT MSVC)
    ADD_CXX_OPTIONS("-Wall")
  ENDIF(NOT MSVC)
  IF(NOT WITH_OCC)
    ADD_CXX_OPTIONS("-Werror=overloaded-virtual")
  ENDIF()
  ADD_CXX_OPTIONS("-Wextra -Wno-unused-parameter")
  ADD_CXX_OPTIONS("-Werror=implicit-function-declaration")
  # should be supported only by Clang. The last statement is important, otherwise nothing compiles ...
  ADD_CXX_OPTIONS("-Werror=conversion -Wno-sign-conversion -Wno-error=sign-conversion Wno-shorten-64-to-32 -Wno-error=shorten-64-to-32")
  # ADD_C_OPTIONS("-Wno-error=shorten-64-to-32") # for clang
  
  IF((NOT WITH_MECHANISMS) AND (NOT SWIG_PROJECT))
    ADD_CXX_OPTIONS("-Werror=missing-declarations")
  ENDIF()
  ADD_CXX_OPTIONS("-Werror=switch-bool")
  ADD_CXX_OPTIONS("-Werror=logical-not-parentheses")
  ADD_CXX_OPTIONS("-Werror=sizeof-array-argument")
  ADD_CXX_OPTIONS("-Werror=bool-compare")
  ADD_CXX_OPTIONS("-Werror=array-bounds")
  ADD_CXX_OPTIONS("-Werror=format-invalid-specifier")
  ADD_CXX_OPTIONS("-Werror=type-limits")
  
  ADD_CXX_OPTIONS("-Wodr")
  
  IF(NOT CXXVERSION STRLESS "201102L" AND DEV_MODE)
    ADD_CXX_OPTIONS("-Wsuggest-final-types")
    ADD_CXX_OPTIONS("-Wsuggest-final-methods")
    
    IF(NOT SWIG_PROJECT)
      ADD_CXX_OPTIONS("-Wzero-as-null-pointer-constant")
    ENDIF()
  ENDIF()
  
  # Compiler Specific
  ADD_CXX_OPTIONS("-diag-disable 654" "Intel")
  IF(NOT ICCOK)
    ADD_CXX_OPTIONS("-D__aligned__=ignored" "Intel")
  ENDIF(NOT ICCOK)
  
  ADD_CXX_OPTIONS("-Wno-string-plus-int" "Clang")
  ADD_CXX_OPTIONS("-Werror=unreachable-code" "Clang")

 IF(USE_SANITIZER MATCHES "asan")
   APPEND_CXX_FLAGS("-fsanitize=leak -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
 ELSEIF(USE_SANITIZER MATCHES "msan")
   APPEND_CXX_FLAGS("-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer")
 ENDIF(USE_SANITIZER MATCHES "asan")

 IF(USE_LIBCXX)
  APPEND_CXX_FLAGS("-stdlib=libc++ -I${USE_LIBCXX}/include -I${USE_LIBCXX}/include/c++/v1")
  SET(_LIBCXX_FLAGS_TO_ADD "-L${USE_LIBCXX}/lib -lc++abi -Wl,-rpath,${USE_LIBCXX}/lib")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${_LIBCXX_FLAGS_TO_ADD}")
 ENDIF(USE_LIBCXX)

ENDIF(CMAKE_CXX_COMPILER_WORKS)
