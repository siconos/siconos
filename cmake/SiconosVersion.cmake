# --- set siconos current version ---
# This file is also required for examples.
set(MAJOR_VERSION 4)
set(MINOR_VERSION 2)
set(PATCH_VERSION 0)
set(SICONOS_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION}")

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

set(SO_current 5)
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
