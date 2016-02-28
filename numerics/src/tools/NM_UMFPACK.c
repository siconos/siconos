/* Siconos-Numerics, Copyright INRIA 2005-2016
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

#include "NumericsMatrix_private.h"
#include "NumericsMatrix.h"
#include "NumericsOptions.h"

#include "debug.h"
#include "misc.h"

#ifdef WITH_UMFPACK

NM_UMFPACK_WS* NM_UMFPACK_factorize(NumericsMatrix* A)
{
  NumericsSparseLinearSolverParams* params = NM_linearSolverParams(A);

  if (params->solver_data)
  {
    return (NM_UMFPACK_WS*) params->solver_data;
  }

  params->solver_data = (NM_UMFPACK_WS*)calloc(1, sizeof(NM_UMFPACK_WS));
  NM_UMFPACK_WS* umfpack_ws = (NM_UMFPACK_WS*) params->solver_data;

  UMFPACK_FN(defaults) (umfpack_ws->control);

  umfpack_ws->control[UMFPACK_PRL] = verbose;
/* TODO UMFPACK_PIVOT_TOLERANCE, UMFPACK_ORDERING, UMFPACK_SCALE
 * UMFPACK_DROPTOL, UMFPACK_STRATEGY, UMFPACK_IRSTEP*/

  CSparseMatrix* C = NM_csc(A);

  csi status;

  status = UMFPACK_FN(symbolic) (C->m, C->n, C->p, C->i, C->x, &(umfpack_ws->symbolic), umfpack_ws->control, umfpack_ws->info);

  if (status)
  {
    umfpack_ws->control[UMFPACK_PRL] = 1;
    UMFPACK_FN(report_status) (umfpack_ws->control, status);
    return NULL;
  }

  status = UMFPACK_FN(numeric) (C->p, C->i, C->x, umfpack_ws->symbolic, &(umfpack_ws->numeric), umfpack_ws->control, umfpack_ws->info);

  if (status)
  {
    umfpack_ws->control[UMFPACK_PRL] = 1;
    UMFPACK_FN(report_status) (umfpack_ws->control, status);
    return NULL;
  }

  umfpack_ws->wi = (csi*)malloc(C->n * sizeof(csi));

  csi size_wd;
  if (umfpack_ws->control[UMFPACK_IRSTEP] > 0)
  {
    size_wd = 5 * C->n;
  }
  else
  {
    size_wd = C->n;
  }
  umfpack_ws->wd = (double*)malloc(size_wd * sizeof(double));

  umfpack_ws->x = (double*)malloc(C->n * sizeof(double));

  return umfpack_ws;
}



void NM_UMFPACK_free(void* p)
{
  assert(p);
  NumericsSparseLinearSolverParams* params = (NumericsSparseLinearSolverParams*) p;
  assert(params);
  NM_UMFPACK_WS* umfpack_ws = (NM_UMFPACK_WS*) params->solver_data;
  assert(umfpack_ws);

  UMFPACK_FN(free_symbolic) (&(umfpack_ws->symbolic));
  UMFPACK_FN(free_numeric) (&(umfpack_ws->numeric));

  if (umfpack_ws->wi)
  {
    free(umfpack_ws->wi);
    umfpack_ws->wi = NULL;
  }

  if (umfpack_ws->wd)
  {
    free(umfpack_ws->wd);
    umfpack_ws->wd = NULL;
  }

  if (umfpack_ws->x)
  {
    free(umfpack_ws->x);
    umfpack_ws->x = NULL;
  }

  /* Here we free umfpack_ws ...  */
  free(umfpack_ws);
  params->solver_data = NULL;

}

void NM_UMFPACK_extra_display(NM_UMFPACK_WS* umfpack_ws)
{
  if (verbose > 2)
  {
    UMFPACK_FN(report_info) (umfpack_ws->control, umfpack_ws->info);

    if (verbose > 3)
    {
      UMFPACK_FN(report_control) (umfpack_ws->control);
    }
  }
  else if (verbose > 1)
  {
    if (umfpack_ws->control[UMFPACK_IRSTEP] > 0)
    {
      printf("UMFPACK : backward error estimate omega1 %g\n", umfpack_ws->info[UMFPACK_OMEGA1]);
      printf("UMFPACK : backward error estimate omega2 %g\n", umfpack_ws->info[UMFPACK_OMEGA2]);
    }
    printf("UMFPACK : solve FLOPS %g\n", umfpack_ws->info[UMFPACK_SOLVE_FLOPS]);
    printf("UMFPACK : solve time %g\n", umfpack_ws->info[UMFPACK_SOLVE_TIME]);
    printf("UMFPACK : wall time %g\n", umfpack_ws->info[UMFPACK_SOLVE_WALLTIME]);

  }
}

#endif /* WITH_UMFPACK */
