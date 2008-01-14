/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
#include <math.h>
#include "NSSTools.h"
#include "LA.h"

int filter_result_LCP_block(SparseBlockStructuredMatrix *blmat, double *q , double *z , double tol, int chat, double *w)
{
  double error, normq;
  double a1, b1;
  int i, j, incx, incy;
  int n, nbbl, nbblrow;
  int rowprecbl, rowcurbl, indicrow, rowsize;
  int colprecbl, colcurbl, indiccol, colsize;

  double *adrcurbl;

  nbbl = blmat->nbblocks;
  nbblrow = blmat->size;
  n = 0;
  for (i = 0 ; i < nbblrow ; i++) n += blmat->blocksize[i];

  a1 = 1.;
  b1 = 1.;

  incx = 1;
  incy = 1;
  DCOPY(n , q , incx , w , incy);

  rowprecbl = -1;
  rowsize = 0;
  indicrow = 0;

  for (i = 0 ; i < nbbl ; i++)
  {
    rowcurbl = blmat->RowIndex[i];
    colcurbl = blmat->ColumnIndex[i];
    if (rowcurbl != rowprecbl)   /*  new row  */
    {
      indiccol = 0;
      for (j = 0 ; j < colcurbl ; j++) indiccol += blmat->blocksize[j];
      colprecbl = colcurbl;

      indicrow += rowsize;
      for (j = rowprecbl + 1 ; j < rowcurbl ; j++) indicrow += blmat->blocksize[j];

      rowprecbl = rowcurbl;
      rowsize = blmat->blocksize[rowcurbl];
    }
    for (j = colprecbl ; j < colcurbl ; j++) indiccol += blmat->blocksize[j];
    colprecbl = colcurbl;

    rowsize = blmat->blocksize[rowcurbl];
    colsize = blmat->blocksize[colcurbl];
    adrcurbl = blmat->block[i];

    DGEMV(LA_NOTRANS , rowsize , colsize , a1 , adrcurbl , rowsize , &z[indiccol] ,
          incx , b1 , &w[indicrow] , incy);
  }

  error = 0.;
  for (i = 0 ; i < n ; i++)
  {
    if (z[i] < 0.0)
    {
      error += -z[i];
      if (w[i] < 0.0) error += z[i] * w[i];
    }
    if (w[i] < 0.0) error += -w[i];
    if ((z[i] > 0.0) && (w[i] > 0.0)) error += z[i] * w[i];
  }

  incx  = 1;
  normq = DNRM2(n , q , incx);

  error = error / normq;
  if (error > tol)
  {
    if (chat > 0) printf(" Wrong LCP block result , error = %g \n", error);
    return 1;
  }
  else return 0;

}
