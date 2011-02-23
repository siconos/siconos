/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef NUMERICSGENERICMECHANICALPROBLEM_H
#define NUMERICSGENERICMECHANICALPROBLEM_H

#include "NumericsMatrix.h"
#include "SolverOptions.h"
/* void * solverFC3D; */
/* void * solverEquality; */
/* void * solverLCP; */
/* void * solverMLCP; */

/** GenericMechanicalProblem
    \param numberOfBlockLine The number of line of blocks.
    \param M A sparse blocks matrix.
    \param q A dense vector.
 */

/** Remark:

    The M and q contains de matrices of the problem. The sub problems (problems) has also a M and q member usfull for the computation of the local error.


 */
#define GMP_MAX_SIZE_LOCAL 30

typedef struct _listNumericsProblem
{
  int type;
  void * problem;
  double *q;/*a pointer on the q of the problem*/
  int size;/*size of the local problem.(needed because of dense case)*/
  struct _listNumericsProblem * nextProblem;
  struct _listNumericsProblem * prevProblem;
} listNumericsProblem;


/*ONLY q and M must be allocated/free by the users, the others fields are private:
  DO NOT FILL THIS STRUCTURE BY YOURSELF, BUT USE THE buildEmptyGenericMechanicalProblem, addProblem and freeGenericMechanicalProblem FUNCTIONS.
 */
typedef struct
{
  /*Number of line of blocks.*/
  /*PRIVATE: manage by addProblem.*/
  int size;
  /*must be set by the user.*/
  NumericsMatrix* M;
  /*must be set by the user.*/
  double* q;
  /*PRIVATE: manage by addProblem.*/
  listNumericsProblem *firstListElem;
  /*PRIVATE: manage by addProblem.*/
  listNumericsProblem *lastListElem;
  //  void * * problems;
} GenericMechanicalProblem;




#endif
