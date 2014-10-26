# apply misc flags

IF(CMAKE_C_COMPILER)

 IF(C_HAVE_C99 AND NOT MSVC)
  APPEND_C_FLAGS("-std=c99")
 ELSEIF(C_HAVE_C99 AND NOT MSVC)
  IF(C_HAVE_XC99 AND NOT MSVC)
   APPEND_C_FLAGS("-xc99")
  ENDIF(C_HAVE_XC99 AND NOT MSVC)
 ENDIF(C_HAVE_C99 AND NOT MSVC)

 IF(C_HAVE_WALL)
  APPEND_C_FLAGS("-Wall")
 ENDIF(C_HAVE_WALL)

 IF(C_HAVE_WEXTRA)
  APPEND_C_FLAGS("-Wextra -Wno-unused-parameter")
 ENDIF(C_HAVE_WEXTRA)

 IF(C_HAVE_IMPL)
  APPEND_C_FLAGS("-Werror=implicit-function-declaration")
 ENDIF(C_HAVE_IMPL)

 IF(C_HAVE_UNREACH AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
  APPEND_C_FLAGS("-Werror=unreachable-code")
 ENDIF(C_HAVE_UNREACH AND NOT "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")

 IF(C_HAVE_CONV)
  APPEND_C_FLAGS("-Werror=conversion -Wno-sign-conversion")
  APPEND_C_FLAGS("-Wno-error=sign-conversion -Wno-error=shorten-64-to-32")
 ENDIF(C_HAVE_CONV)

 # too many errors right now ...
 #IF(C_HAVE_MISS)
 #  APPEND_C_FLAGS("-Wmissing-prototypes")
 #ENDIF(C_HAVE_MISS)

ENDIF(CMAKE_C_COMPILER)

IF(CMAKE_CXX_COMPILER)

 IF(NOT CMAKE_COMPILER_IS_GNUCXX)
  IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
   # Disable warnings with intel compiler due (mainly) to visitors visit function overloading
   IF(CXX_HAVE_DIAG_DISABLE_654)
    APPEND_CXX_FLAGS("-diag-disable 654")
   ENDIF(CXX_HAVE_DIAG_DISABLE_654)
   # Error on intel compiler, see: http://software.intel.com/en-us/forums/showthread.php?t=65041
   # This issue have been solved with ICC >= 12.1
   if(NOT ICCOK)
    if(CXX_HAVE_D__ALIGNED__IGNORED)
     APPEND_CXX_FLAGS("-D__aligned__=ignored")
    endif(CXX_HAVE_D__ALIGNED__IGNORED)
   endif(NOT ICCOK)
  endif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
 endif(NOT CMAKE_COMPILER_IS_GNUCXX)
 # way too verbose with MSVC
 IF(CXX_HAVE_WALL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Wall")
 ENDIF(CXX_HAVE_WALL AND NOT MSVC)
 IF(CXX_HAVE_WEXTRA)
  APPEND_CXX_FLAGS("-Wextra -Wno-unused-parameter")
 ENDIF(CXX_HAVE_WEXTRA)
 IF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   # stupid warning
  APPEND_CXX_FLAGS("-Wno-string-plus-int")
 ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")

 IF(CXX_HAVE_WOVERLOADED_VIRTUAL)
   APPEND_CXX_FLAGS("-Woverloaded-virtual")
 ENDIF(CXX_HAVE_WOVERLOADED_VIRTUAL)

 IF(CXX_HAVE_IMPL AND NOT MSVC)
  APPEND_CXX_FLAGS("-Werror=implicit-function-declaration")
 ENDIF(CXX_HAVE_IMPL AND NOT MSVC)

 IF(CXX_HAVE_UNREACH AND NOT MSVC AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    APPEND_CXX_FLAGS("-Werror=unreachable-code")
 ENDIF(CXX_HAVE_UNREACH AND NOT MSVC AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

 IF(CXX_HAVE_CONV AND NOT MSVC)
   APPEND_CXX_FLAGS("-Werror=conversion -Wno-sign-conversion")
   APPEND_CXX_FLAGS("-Wno-error=sign-conversion -Wno-error=shorten-64-to-32")
 ENDIF(CXX_HAVE_CONV AND NOT MSVC)

 IF(CXX_HAVE_MISS AND NOT MSVC)
    APPEND_CXX_FLAGS("-Wmissing-declarations")
 ENDIF(CXX_HAVE_MISS AND NOT MSVC)

ENDIF(CMAKE_CXX_COMPILER)


