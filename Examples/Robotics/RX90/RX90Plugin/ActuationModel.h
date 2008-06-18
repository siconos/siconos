
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
** @file ActuationModel/ActuationModel.h
** @author: Pierre-Brice Wieber
** Affiliation(s): INRIA, team BIPOP
** Email(s): Pierre-Brice.Wieber@inria.fr
**
** @brief Declarations of the functions defining the model of actuators
**
**
*/
#ifndef __ActuationModel_h
#define __ActuationModel_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/**
 * Initialization of actuators dynamics for a position q, a velocity qdot and a time of sample given
 * @param[in] t (double) time
 * @param[in] q (double vector, size = NDOF) position
 * @param[in] qdot (double vector, size = NDOF) velocity
 * @param[in] contactState (int vector, size = NCONT) state of contacts
 * @param[in] NDOF (int)
 * @param[in] NCONT (int)
 *
 * @param[out] state (double vector) state of actuators switchs
 * @param[out] z (double vector)
 */
extern
#ifdef __cplusplus
"C"
#endif
void actuationInitialisation(double *t, double *q, double *qdot, int *contactState, int *NDOF, int *NCONT, double *state, double *z);

/**
 * Compute the State of actuators switchs
 * @param[in] t (double) time
 * @param[in] x (double vector, size = 2*NDOF + 2*Nmuscles + 2) system state \f$x = (q, \dot{q}, z)\f$
 * @param[in] oldState (double vector) state of actuators switchs
 * @param[in] NDOF (int) number of degrees of freedom
 * @param[in] x_size (int) size of x
 * @param[in] state_size (int) size of oldState
 *
 * @param[out] newState (double vector, size = nAct) new state of actuators switchs
 */
extern
#ifdef __cplusplus
"C"
#endif
void actuationStateReset(double *t, double *x, double *oldState, int *NDOF, int *x_size, int *state_size, double *newState);

/**
 * Compute the state of the actuators switchs
 * @param[in] t (double) time
 * @param[in] x (double vector, size = 2*NDOF + ??) system state \f$x = (q, \dot{q}, z)\f$
 * @param[in] actuationState (double vector)
 * @param[in] NDOF (int) number of degrees of freedom
 * @param[in] x_size (int) size of vector x
 * @param[in] state_size (int) size of vector actuationState
 *
 * @param[out] events (double vector)
 * @param[out] events_size (int) size of the events vector
 */
extern
#ifdef __cplusplus
"C"
#endif
void actuationEventDetection(double *t, double *x, double *actuationState, int *NDOF, int *x_size, int *state_size, double *events, int *events_size);

/**
 * Handling of Actuation Events
 * @param[in] t (double) time
 * @param[in] x (double vector, size = 2*NDOF ) system state \f$x = (q, \dot{q}, z)\f$
 * @param[in] events (int vector) state of the actuators switchs
 * @param[in] oldState (double vector) state of actuators switchs
 * @param[in] NDOF (int) number of degrees of freedom
 * @param[in] x_size (int) size of the x vector
 * @param[in] state_size (int) size of the oldState vector
 * @param[in] NCONT (int) size of the contactState vector
 * @param[in] contactState (int vector, size = NCONT) state of contacts
 *
 * @param[out] newx (double vector, size = 2*NDOF ) new system state \f$x = (q, \dot{q}, z)\f$
 * @param[out] newState (int vector) new state of actuators switchs
 * @param[out] warmStart (boolean)
 */
extern
#ifdef __cplusplus
"C"
#endif
void actuationEventHandling(double *t, double *x, int *events, double *oldState, int *NDOF, int *x_size, int *state_size, int *NCONT, int *contactState, double *newx, double *newState, int *warmStart);

/**
 * Compute the generalized Torques
 * @param[in] t (double) time
 * @param[in] q (double vector, size = NDOF) position
 * @param[in] qdot (double vector, size = NDOF) velocity
 * @param[in] z (double vector)
 * @param[in] state (double vector)
 * @param[in] NDOF (int)
 * @param[in] NCONT (int)
 * @param[in] z_size (int)
 *
 * @param[out] zdot (double vector)
 * @param[out] torques (double vector, size = NDOF) generalized torques
 */
extern
#ifdef __cplusplus
"C"
#endif
void actuationDynamics(double *t, double *q, double *qdot, double *z, double *state, int *NDOF, int *NCONT, int *z_size, double *zdot, double *torques);
#endif /* __ActuationModel_h */
