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
** @file ActuationModel/NoDynamics/TaskFunctionControl/Trajectory.h
** @author: Florence Billet
** Affiliation(s): INRIA, team BIPOP
** Email(s): Florence.Billet@inria.fr
**
** @brief Header for control law and trajectory in NoDynamics
**
*/
#ifndef __TaskFunctionDefinitoin_trajectory_h
#define __TaskFunctionDefinitoin_trajectory_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * Compute the position, velocity and acceleration desired at a given time t
 *  @param[in] t (double) time
 *
 *  @param[out] position (double vector, size = NDOF) desired position in task space
 *  @param[out] velocity (double vector, size = NDOF) desired velocity in task space
 *  @param[out] acceleration (double vector, size = NDOF) desired acceleration
 *  in task space
 *  @param[out]contacts (int vector, variable size <= NCONTACTS ) number of the
 *  lines of contact vector which correspond to effectives contacts at t.
 */
void trajectory(double * t, double * q, double * qdot, double * qddot, int * contacts);

#endif
