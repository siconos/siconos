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
*/
/**
** @file LagrangianModel/Bip/LagrangianModel.h
** @author: Pierre-Brice Wieber
** Affiliation(s): INRIA, team BIPOP
** Email(s): Pierre-Brice.Wieber@inria.fr
**
** @brief Compute the complete Dynamics
**
**
*/
#ifndef __LagrangianModel_h
#define __LagrangianModel_h

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif

/* A mettre dans la generation maple !!! */
#define NDOF 7
#define NTAGS (7+1)
#define NCONT 2

/**
 * Compute the matrix of contact for a given bipede state
 *
 * @param[out] CC Matrix of contact <em> dim NCONTx3(xyz) </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void
Contact(double *CC, double *q);

/**
 * Compute the matrix of contact jacobian for a given bipede state
 *
 * @param[out] CJ Matrix of Contact Jacobian <em> dim (3xNCONT)*NDOF </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void
ContactJacobian(double *CJ, double *q);

/**
 * Compute the vector of contact hessians for a given bipede state
 *
 * @param[out] H Vector of Contact Hessians <em> dim NCONT*3 </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 * @param[in] qdot Articular Velocity State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void
ContactHessian(double *H, double *q, double *qdot);

/**
 * Compute the matrix of Non Linear Effect (Coriolis + Gravity)
 * for a given bipede state
 *
 * @param[out] N Matrix of Coriolis <em> dim NDOF </em>
 * @param[in] q joint State Vector <em> dim NDOF </em>
 * @param[in] qdot Articular Velocity State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void
NLEffects(double *N, double *q, double *qdot);

/**
 * Compute the matrix of inertia for a given bipede state
 *
 * @param[out] M Matrix of Inertia <em> dim NDOFxNDOF </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */


SICONOS_EXPORT void JacobianNLEffects(double *NJ, double *q, double *qdot);

SICONOS_EXPORT void JacobianVelocityNLEffects(double *NVJ, double *q, double *qdot);




SICONOS_EXPORT void
Inertia(double *M, double *q);

/**
 * Compute the matrix of tags for a given bipede state
 * This matrix is made up of caracteristic points of bipede
 * in the absolute referential. The end of matrix T contains
 * the coordinate of the biped mass center
 *
 * @param[out] T Matrix of contact <em> dim 60 = NTAGSx3(xyz) </em>
 * @param[in] q Joint State Vector <em> dim NDOF </em>
 */
SICONOS_EXPORT void
Tags(double *T, double *q);



#endif /* __LagrangianModel_h */
