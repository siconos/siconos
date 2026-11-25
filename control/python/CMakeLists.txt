# SPDX-FileCopyrightText: INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

#================================================================
# cmake utilities to build python-bindings for control component
#================================================================

## PB11 wrap for control NOT YET IMPLEMENTED - TO BE DONE


if(WITH_TESTING)
  if(NOT HAS_FORTRAN)
    list(APPEND python_excluded_tests tests/test_smc.py)
  endif()
  # NOT YET IMPLEMENTED
  # build_python_tests(
  #   # DEPS ${COMPONENT} # some plugins in tests need to be linked with kernel
  #   EXCLUDE ${python_excluded_tests})
endif()
