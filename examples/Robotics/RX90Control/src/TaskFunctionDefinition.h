/* Copyright (C) INRIA 1999-2008
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
*/
/**
** @file ActuationModel/NoDynamics/TaskFunctionControl/TaskFunctionDefinition.h
** @author: Florence Billet
** Affiliation(s): INRIA, team BIPOP
** Email(s): Florence.Billet@inria.fr
**
** @brief Task Function implementation
**
*/
#ifndef __taskfunctiondefinition_h
#define __taskfunctiondefinition_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif



/**
 * Compute the vector s(q) matching to the task function for a given bipede state q
 *
 * @param[out] T vector matching to the output function <em> dim NDOF </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void TaskFunction(double *T, double *q);

/**
 * Compute the jacobian matrix of the task function for a given bipede state q
 *
 * @param[out] J Jacobian matrix of the output function <em> dim NDOFxNDOF </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void TaskJacobian(double J[441], double q[21]);

/**
 * Compute the matrix of Non Linear Effect (Coriolis + Gravity)
 * for a given bipede state
 *
 * @param[out] H Matrix of Coriolis <em> dim NDOF </em>
 * @param[in] q joint State Vector <em> dim NDOF </em>
 * @param[in] qdot Articular Velocity State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void TaskNLEffects(double H[21], double q[21], double qdot[21]);
#endif

