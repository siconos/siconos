/* Copyright (C) INRIA 1999-2005
**
** This program is free software; you can redistribute it and/or modify it
** under the terms of the GNU General Public License version 2 as published
** by the Free Software Foundation.
**
** This program is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
** Public License for more details.
**
** You should have received a copy of the GNU General Public License along
** with this program; if not, write to the Free Software Foundation, Inc.,
** 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**
** Author(s): Pierre-Brice Wieber
** Affiliation(s): INRIA, team BIPOP
** Email(s): Pierre-Brice.Wieber@inria.fr
**
** Description:
**
** Modifications:
** $Log: LagrangianModel.h,v $
** Revision 3.0.0.1  2005/02/08 13:05:22  rpissard
** version start HuMAnS
**
*/
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>

SICONOS_EXPORT void NLEffects(double *N, double* q, double *qdot);

SICONOS_EXPORT void Inertia(double *M, double *q);

SICONOS_EXPORT void JacobianQNLEffects(double *jaco, double* q, double *qdot);

SICONOS_EXPORT void JacobianVNLEffects(double *jaco, double* q, double *qdot);

SICONOS_EXPORT void Contact(double *CC, double* q);

SICONOS_EXPORT void ContactJacobian(double *CJ, double *q);
