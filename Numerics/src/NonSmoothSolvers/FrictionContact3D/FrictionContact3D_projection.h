/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#ifndef FRICTIONCONTACT3DProjection_H
#define FRICTIONCONTACT3Projection_H

/*!\file FrictionContact3D_projection.h
  \brief Typedef and functions declarations related to projection solver for 3 dimension frictional contact problems

  Each solver must have 4 functions in its interface:
  - initialize: link global static variables to the considered problem (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

  We consider a "global" (ie for several contacts) problem, used to initialize the static global variables.
  Then a "local" (ie for one contact => size = 3) problem is built (update function) and solved (solve function).

  Two different storages are available for M: dense and sparse block.

  \author INRIA Siconos Team

*/
#include "NumericsMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** Initialize friction-contact 3D projection
      \param dim. of the global problem
      \param matrix M of the global problem
      \param vector q of the global problem
      \param vector of the friction coefficients
  */
  void frictionContact3D_projection_initialize(int, const NumericsMatrix*const, const double*const, const double*const);

  /** Initialize friction-contact 3D projection with regularization
      \param dim. of the global problem
      \param matrix M of the global problem
      \param vector q of the global problem
      \param vector of the friction coefficients
  */
  void frictionContact3D_projection_initialize_with_regularization(int, const NumericsMatrix*const, const double*const, const double*const);

  /** Update friction-contact 3D projection solver: formalize local problem for one contact.
      \param number (position in global matrix) of the considered contact
      \param global reaction (only the block corresponding to the current contact will be modified,
      the rest is used to formalize the local problem)
  */
  void frictionContact3D_projection_update(int, double*);

  /** solve friction-contact 3D problem with projection assuming that M is diagonal
      \param number (position in global matrix) of the considered contact
      \param dim. of the global problem
      \param global reaction (only the block corresponding to the current contact will be modified,
      \param vector of int parameters (max iteration numnber ...)
      \param vector of double parameters (tolerance ...)
  */
  void frictionContact3D_projectionWithDiagonalization_solve(int, int, double*, int*, double*);

  /** solve friction-contact 3D problem with projection on the Cone
      \param number (position in global matrix) of the considered contact
      \param dim. of the global problem
      \param global reaction (only the block corresponding to the current contact will be modified,
      \param vector of int parameters (max iteration numnber ...)
      \param vector of double parameters (tolerance ...)
  */
  void frictionContact3D_projectionOnCone_solve(int, int, double*, int*, double*);

  /** solve friction-contact 3D problem with projection on the Cone and regularization

     \param number (position in global matrix) of the considered contact
     \param dim. of the global problem
     \param global reaction (only the block corresponding to the current contact will be modified,
     \param vector of int parameters (max iteration numnber ...)
     \param vector of double parameters (tolerance ...)
  */
  void frictionContact3D_projectionOnCone_with_regularization_solve(int, int, double*, int*, double*);

  /** solve friction-contact 3D problem with projection on the Cone with local
      iteration up to convergence of the local problem
      \param number (position in global matrix) of the considered contact
      \param dim. of the global problem
      \param global reaction (only the block corresponding to the current contact will be modified,
      \param vector of int parameters (max iteration numnber ...)
      \param vector of double parameters (tolerance ...)
  */
  void frictionContact3D_projectionOnConeWithLocalIteration_solve(int, int, double*, int*, double*);

  /** free memory for friction contact 3D projection solver */
  void frictionContact3D_projection_free();


  /** Implement the projection on the Coulomb cone
      \param a pointrt to the 3D vector of Reaction
      \param the coefficient of friction
  */
  void projectionOnCone(double* reaction, double  mu);



#ifdef __cplusplus
}
#endif

#endif
