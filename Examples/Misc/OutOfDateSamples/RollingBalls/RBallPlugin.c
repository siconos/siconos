/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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
#include <stdio.h>

const double R = 0.1; // ball radius
const double m = 1; // ball mass
const double g = 9.8; // gravity
const double h = -0.1; // altitude of the plan

extern "C" void vectorField(int *sizeOfX, double *time, double *x, double *xdot)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output patarmeter : xdot (pointer to Xdot vector)
   */

  printf("Call of the function 'vectorField' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void ballMass(int *sizeOfq, double *time, double *q, double* mass)
{
  /* input parameter : sizeOfq (size of the vector q); time ; q (pointer to q vector);
   * output patarmeter : mass (pointer to mass matrix)
   */
  //  printf("Call of the function 'ballMass' of the basic plugin.\nYou have to implement this function.\n");

  int n = *sizeOfq;

  double I = (3 / 5.0) * m * R * R;

  int i;

  // initialisation of the Mass matrix
  for (i = 0; i < n ^ 2; i++)
  {
    mass[i] = 0.0;
  }

  mass[0] = m;
  mass[4] = m;
  mass[8] = I;
}

extern "C" void ballQNLInertia(int *sizeOfq, double *q, double *velocity, double *Q)
{
  /* input parameter : sizeOfq (size of the vector q); q (pointer to q vector); velocity (pointer to velocity vector);
   * output patarmeter : Q (pointer to Q vector)
   */

  //  printf("Call of the function 'ballQNLInertia' of the basic plugin.\nYou have to implement this function.\n");
  int i;
  int n = *sizeOfq;

  for (i = 0; i < n; i++)
  {
    Q[i] = 0.0;
  }
}

extern "C" void ballFInt(int *sizeOfq, double *time, double *q, double *velocity, double *fInt)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output patarmeter : fInt (pointer to Fint *vector)
   */

  //  printf("Call of the function 'ballFInt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;

  for (i = 0; i < n; i++)
  {
    fInt[i] = 0.0;
  }
}

extern "C" double FextFunction(double time)
{
  double t = time;

  double res = -0.0;
  return res;
}


extern "C" void ballFExt(int *sizeOfq, double *time, double *fExt)
{
  /* input parameter : sizeOfq (size of the vector q); time;
   * output patarmeter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'ballFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  double t = *time;

  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }

  //  fExt[0] = -m*g + FextFunction(t);
}

extern "C" void groundFExt(int *sizeOfq, double *time, double *fExt)
{
  /* input parameter : sizeOfq (size of the vector q); time;
   * output patarmeter : fExt (pointer to Fext vector)
   */

  //  printf("Call of the function 'ballFExt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  double t = *time;

  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }
}

extern "C" void ballJacobianX(double *time, double *x, double *jacob)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output patarmeter : jacob (pointer to JacobianX matrix)
   */

  //printf("Call of the function 'ballJacobianX' of the basic plugin.\nYou have to implement this function.\n");
}


extern "C" void ballJacobianQFInt(int *sizeOfq, double *time, double *q, double *velocity, double *jacob)
{
  /* input parameter : sizeOfq (size of the vector q); time; q (pointer to q vector); velocity (pointer to velocity vector);
   * output patarmeter : jacob (pointer to JacobianCoordFint *matrix)
   */

  //  printf("Call of the function 'ballJacobianCoordFInt' of the basic plugin.\nYou have to implement this function.\n");
  int i;
  int n = *sizeOfq;
  for (i = 0; i < n; i++)
  {
    jacob[i] = 0.0;
  }
}

extern "C" void ballJacobianVelocityFInt(int *sizeOfq, double *time, double *q, double *velocity, double *jacob)
{
  //  printf("Call of the function 'ballJacobianVelocityFInt' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  for (i = 0; i < n; i++)
  {
    jacob[i] = 0.0;
  }
}

extern "C" void ballJacobianQQNLInertia(int *sizeOfq, double *q, double *velocity, double *jacob)
{
  //printf("Call of the function 'ballJacobianCoordQ' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  for (i = 0; i < n; i++)
  {
    jacob[i] = 0.0;
  }
}

extern "C" void ballJacobianVelocityQNLInertia(int *sizeOfq, double *q, double *velocity, double *jacob)
{
  //  printf("Call of the function 'ballJacobianVelocityQ' of the basic plugin.\nYou have to implement this function.\n");

  int i;
  int n = *sizeOfq;
  for (i = 0; i < n; i++)
  {
    jacob[i] = 0.0;
  }
}

extern "C" void ballOutput(int *sizeOfX, double* x, double *time, int *sizeOfLambda, double* lambda, double* y)
{
  /* input parameter : sizeOfX (size of the vector X);
   * x (pointer to x vector);
   * time;
   * sizeOfLambda(size of the vector Lambda);
   * lambda (pointer to lambda vector);
   * output parameter : y (pointer to vector y )
   */
  //  printf("Call of the function 'ballOutput' of the basic plugin.\nYou have to implement this function.\n");

  y[1] = x[1] - (2 * R) - x[4]; // - h;
  printf("$$ extern \"C\" void ballOutput\n");
}


extern "C" void ballInput(int *sizeOfX, double* x, double *time, int *sizeOfLambda, double* lambda, double* r)
{
  /* input parameter : sizeOfX (size of the vector X); x (pointer to x vector); time; sizeOfLambda(size of the vector Lambde); lambda (pointer to lambda vector)
   * output patarmeter : r (pointer to vector r )
   */
  //  printf("Call of the function 'ballInput' of the basic plugin.\nYou have to implement this function.\n");

  r[0] = lambda[0];
  r[1] = 0;
  r[2] = 0;
  r[3] = -lambda[0];
}

extern "C" void ballU(/* to be defined */)
{
  printf("Call of the function 'ballU' of the basic plugin.\nYou have to implement this function.\n");
}

extern "C" void ballF(/* to be defined */)
{
  printf("Call of the function 'ballF' of the basic plugin.\nYou have to implement this function.\n");
}

