// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file ActuationModel/NoDynamics/ActuationDynamics.cpp
// @author Pierre-Brice Wieber, Christine Azevedo
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief Compute the generalized Torques
//

// Description:
//
// Modifications:
//
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <iostream>

extern "C" {
#include "NoDynamics.h"
}

SICONOS_EXPORT void actuationDynamics(double *t, double *q, double *qdot, double *z, double *state, int *NDOF, int *NCONT, int *z_size, double *zdot, double *torques)
{
  controlLaw(t, q, qdot, NDOF, NCONT, torques);
  zdot = 0;
}
