# =====================================
# Some macros to check C compiler,
# check/add some flags, get version ...
# =====================================

include(CheckCCompilerFlag)

# --- Append some options (list OPT) to cxx compiler flags. ---
# Update :
# - CMAKE_C_FLAGS
# -
# parameters :
# OPT : the flags
# extra parameters : a list of compilers ids
# (see http://www.cmake.org/cmake/help/v3.3/variable/CMAKE_LANG_COMPILER_ID.html#variable:CMAKE_<LANG>_COMPILER_ID)
# If extra parameter is set, the option will be applied
# only for those compilers, if the option is accepted.
# 
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
endmacro()

