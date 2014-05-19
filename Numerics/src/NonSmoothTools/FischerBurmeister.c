/* Siconos-Numerics, Copyright INRIA 2005-2014
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "SiconosBlas.h"
#include "FischerBurmeister.h"
#include "assert.h"

/* Computation of  Fischer-Burmeister function, phi(z,F(z)) = sqrt(z*z + F(z)*F(z)) - z - F(z) */
void phi_FB(int size, double* z, double* F, double* phiVector)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(phiVector != NULL);

  for (int i = 0 ; i < size ; ++i)
  {
    phiVector[i] = sqrt(z[i] * z[i] + F[i] * F[i]) - (z[i] + F[i]);
  }
}

/* Compute the jacobian of the Fischer function */
void jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* jacobianPhiMatrix)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(jacobianF != NULL);
  assert(jacobianPhiMatrix != NULL);

  /* jacobianPhiMatrix is initialized with jacobianF */
  cblas_dcopy(size * size, jacobianF, 1, jacobianPhiMatrix, 1);

  double ri, ai, bi;

  for (int i = 0; i < size; i++)
  {
    ri = sqrt(z[i] * z[i] +  F[i] * F[i]);
    if (ri > 0.0)
    {
      ai = z[i] / ri - 1.0;
      bi = F[i] / ri - 1.0;
    }
    else
    {
      ai = -1.0;
      bi = -1.0;
    }
    /* jacobianPhiMatrix_ij = delta_ij*ai + bi * jacobianF_ij
       delta_ij being the Kronecker symbol
    */
    cblas_dscal(size, bi, &jacobianPhiMatrix[i], size);
    jacobianPhiMatrix[i * size + i] += ai;
  }
}

/* Computation of the mixed Fischer-Burmeister function */
void phi_Mixed_FB(int sizeEq, int sizeIneq, double* z, double* F, double* phiVector)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(phiVector != NULL);

  int totalSize = sizeEq + sizeIneq;

  for (int i = 0 ; i < sizeEq ; ++i)
  {
    phiVector[i] = F[i];
  }
  for (int i = sizeEq ; i < totalSize ; ++i)
  {
    phiVector[i] =  sqrt(z[i] * z[i] + F[i] * F[i]) - (z[i] + F[i]);
  }
}

/* Compute the jacobian of the mixed Fischer function */
void jacobianPhi_Mixed_FB(int sizeEq, int sizeIneq, double* z, double* F, double* jacobianF, double* jacobianPhiMatrix)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(jacobianF != NULL);
  assert(jacobianPhiMatrix != NULL);

  /* jacobianPhiMatrix is initialized with jacobianF */
  cblas_dcopy((sizeEq + sizeIneq) * (sizeEq + sizeIneq), jacobianF, 1, jacobianPhiMatrix, 1);

  double ri, ai, bi;

  for (int i = sizeEq; i < sizeEq + sizeIneq; i++)
  {
    ri = sqrt(z[i] * z[i] +  F[i] * F[i]);
    if (ri > 0.0)
    {
      ai = z[i] / ri - 1.0;
      bi = F[i] / ri - 1.0;
    }
    else
    {
      ai = -1.0;
      bi = -1.0;
    }
    /* jacobianPhiMatrix_ij = delta_ij*ai + bi * jacobianF_ij
       delta_ij being the Kronecker symbol
    */
    /*        cblas_dscal(size, bi, &jacobianPhiMatrix[i], size);*/

    cblas_dscal(sizeEq + sizeIneq, bi, &jacobianPhiMatrix[i], sizeEq + sizeIneq);
    jacobianPhiMatrix[i * (sizeEq + sizeIneq) + i] += ai;
  }
}

void Jac_F_FB(int n1, int n2, double* z, double* F, double* workV1, double* workV2, double* nabla_F, double* H)
{
  double normi;
  int n = n1 + n2;

  if (n1 > 0)
  {
    //printf("Jac_F_FB: the mixed case needs review and testing -- xhub\n");
    /* H is initialized with nabla_F */
    cblas_dcopy((n2+n1) * (n2+n1), nabla_F, 1, H, 1);
  }
  // constructing the set beta
  // Introduce a tolerance ? -- xhub
  for (int i = n1; i < n; ++i)
  {
    if ((fabs(z[i]) < DBL_EPSILON) && (fabs(F[i]) < DBL_EPSILON))
    {
      workV1[i] = 1.0;
    }
    else
    {
      workV1[i] = 0.0;
    }
  }
  // workV1 = "z" in Facchinei--Pang (2003) p. 808
  // "z_i" = 1 if z_i = w_i = 0.0
  // nabla_F^T.workV1 --> workV2
  cblas_dgemv(CblasColMajor,CblasTrans, n2 , n2 , 1.0 , nabla_F , n2 , &workV1[n1], 1, 0.0 , &workV2[n1], 1);
  for (int i = n1; i < n; ++i)
  {
    if (workV1[i] != 0.0) // i in beta
    {
      normi = sqrt(workV1[i] * workV1[i] + workV2[i] * workV2[i]);
      for (int j = 0; j < n; j++)
      {
        H[j * n + i] = (workV2[i] / normi - 1.0) * nabla_F[j * n + i];
      }
      H[i * n + i] += (workV1[i] / normi - 1.0);

    }
    else // i not in beta
    {
      normi = sqrt(z[i] * z[i] + F[i] * F[i]);
      for (int j = 0; j < n; j++)
      {
        H[j * n + i] = (F[i] / normi - 1.0) * nabla_F[j * n + i];
      }
      H[i * n + i] += (z[i] / normi - 1.0);
    }
  }
}
