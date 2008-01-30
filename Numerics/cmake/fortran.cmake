# this file comes from http://freeeos.sourceforge.net/
#
# cmake/modules/fortran.cmake
#
# Copyright (C) 2006 Alan W. Irwin
#
#
#This file is part of FreeEOS.
#
#This file is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; version 2 of the License.
#
#This file is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Library General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with the file; if not, write to the Free Software
#Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
# Module for determining F77/F95 bindings configuration options

# Options to enable Fortran bindings
option(ENABLE_f77 "Enable f77 bindings" ON)
option(ENABLE_f95 "Enable f95 bindings" ON)

if(ENABLE_f77 OR ENABLE_f95)
  # enable_language does not provide proper support for special compiler
  # flags so specify Fortran as the language in the PROJECT statement
  # instead.
  # enable_language(Fortran)
  # Check for fortran compiler
  if(NOT CMAKE_Fortran_COMPILER)
    message(FATAL_ERROR 
    "fortran compiler not found. Nothing will be built."
    )
  endif(NOT CMAKE_Fortran_COMPILER)

  # Don't compile Fortran 95 binding if compiler doesn't support it
  if(ENABLE_f95 AND NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    message(STATUS "WARNING: " 
    "fortran compiler does not support f90/95. Disabling f95"
    )
    set(ENABLE_f95 OFF CACHE BOOL "Enable f95 bindings" FORCE)
  endif(ENABLE_f95 AND NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

  # Set installation location for f95 modules.
  set(F95_MOD_DIR ${LIB_DIR}/fortran/modules/${PACKAGE})
endif(ENABLE_f77 OR ENABLE_f95)

if(NOT ENABLE_f77 AND NOT ENABLE_f95)
  message(FATAL_ERROR 
  "One of ENABLE_f77 or ENABLE_f95 must be set to ON"
  )
endif(NOT ENABLE_f77 AND NOT ENABLE_f95)
