# --- set siconos current version ---
#
# The release number is written in ONE place: the file VERSION at the root of
# the sources. Calendar versioning YEAR.MONTH[.PATCH] (e.g. 2026.11, then
# 2026.11.1 for a bug-fix release). Use dots, not '-': CMake and Python
# (PEP 440, where 2026-11 means 2026.post11) need a dotted numeric version.
#
# Everything else (libraries, siconos-config-version.cmake, python package,
# scripts, doxygen, sphinx) uses the variables below. codemeta.json (read by
# Software Heritage) is static and is only checked here.
set(_siconos_root "${CMAKE_CURRENT_LIST_DIR}/..")
file(STRINGS "${_siconos_root}/VERSION" _siconos_version LIMIT_COUNT 1)
string(STRIP "${_siconos_version}" _siconos_version)
if(NOT _siconos_version MATCHES "^([0-9]+)\\.([0-9]+)(\\.([0-9]+))?$")
  message(FATAL_ERROR "VERSION file: '${_siconos_version}' is not YEAR.MONTH[.PATCH]")
endif()
set(MAJOR_VERSION ${CMAKE_MATCH_1})
set(MINOR_VERSION ${CMAKE_MATCH_2})
if("${CMAKE_MATCH_4}" STREQUAL "")
  set(PATCH_VERSION 0)
else()
  set(PATCH_VERSION ${CMAKE_MATCH_4})
endif()
set(SICONOS_VERSION "${_siconos_version}")

get_property(_codemeta_checked GLOBAL PROPERTY SICONOS_CODEMETA_CHECKED)
if(EXISTS "${_siconos_root}/codemeta.json" AND NOT _codemeta_checked)
  set_property(GLOBAL PROPERTY SICONOS_CODEMETA_CHECKED TRUE)
  file(READ "${_siconos_root}/codemeta.json" _codemeta)
  string(JSON _codemeta_version ERROR_VARIABLE _codemeta_err GET "${_codemeta}" version)
  if(NOT _codemeta_err STREQUAL "NOTFOUND" OR NOT _codemeta_version STREQUAL SICONOS_VERSION)
    message(WARNING "codemeta.json version '${_codemeta_version}' differs from VERSION "
      "'${SICONOS_VERSION}': update codemeta.json (version, dateModified) before the release.")
  endif()
endif()

### SOVERSION
#
# The SICONOS_VERSION variable above indicates API compatibility,
# whereas the SICONOS_SOVERSION (or SONAME), below, indicates ABI
# compatibility.
#
# To be bumped at each release, by the following rules.  If you are
# not sure, likely API and ABI compatibility have both been
# sacrificed, so simply bump SO_current and set the others to zero.
# If an effort has been made to be backwards compatible on this
# release (e.g. bug fix release), continue with the rules outlined
# below.
#
### RULES for SONAME (borrowed from libtool)
### https://www.gnu.org/software/libtool/manual/html_node/Libtool-versioning.html
#
# If the library source code has changed at all since the last update, then
# increment revision (`c:r:a' becomes `c:r+1:a').
#
# If any interfaces have been added, removed, or changed since the last update,
# increment current, and set revision to 0.
#
# If any interfaces have been added since the last public release, then
# increment age.
#
# If any interfaces have been removed since the last public release, then set
# age to 0.

set(SO_current 8)
set(SO_revision 0)
set(SO_age 0)

# Aggregate variables, to be passed to linker.
# libraries will be named e.g.,
#   libsiconos_kernel.so -> libsiconos_kernel.so.5 -> libsiconos_kernel.so.5.0.0
# Again: this is *not* the software release number!
set(SO_version_info "${SO_current}:${SO_revision}:${SO_age}")
math(EXPR SO_current_minus_age "(${SO_current}) - (${SO_age})")
set(SICONOS_SOVERSION "${SO_current_minus_age}.${SO_revision}.${SO_age}" CACHE STRING "Siconos SONAME")
set(SICONOS_SOVERSION_MAJOR "${SO_current_minus_age}" CACHE STRING "Siconos SONAME current-minus-age")
