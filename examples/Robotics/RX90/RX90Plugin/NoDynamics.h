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
** @file ActuationModel/NoDynamics/NoDynamics.h
** @author: Florence Billet
** Affiliation(s): INRIA, team BIPOP
** Email(s): Florence.Billet@inria.fr
**
** @brief Header for control law and trajectory in NoDynamics
**
*/
#ifndef __noDynamics_h
#define __noDynamics_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * Compute the generalized Torques
 *  @param[in] t (double) time
 *  @param[in] q (double vector, size = NDOF) position
 *  @param[in] qdot (double vector, size = NDOF) velocity
 *  @param[in] NDOF (int) number of degrees of freedom
 *  @param[in] NCONT (int) number of contacts
 *
 *  @param[out] torques (double vector, size = NDOF) generalized Torques
 */
void controlLaw(double * t, double * q, double * qdot, int * NDOF, int * NCONT, double * torques);

#endif
