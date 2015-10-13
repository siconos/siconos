/* Siconos-Numerics, Copyright INRIA 2005-2015
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

#include "SiconosConfig.h"

#ifdef HAVE_PATHFERRIS

#include <limits.h>
#include <stdio.h>
#include <assert.h>

#include "MCP_Interface.h"

#include "Path.h"
#include "PathOptions.h"

#include "Macros.h"
#include "Output_Interface.h"
#include "Options.h"

#include "NonlinearComplementarityProblem.h"
#include "PathAlgebra.h"
#include "SolverOptions.h"

#include "Path_interface.h"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

#if defined(USE_OUTPUT_INTERFACE)
/* callback to register with PATH: output from PATH will go here */
static CB_FUNC(void) messageCB (void *data, int mode, char *buf)
{
  fprintf (stdout, "%s", buf);
} /* messageCB */
#endif

void SN_path_interface(MCP_Interface* restrict mcp_interface, double* restrict z, double* restrict F, int* restrict info , double* restrict dparam, int* restrict iparam)
{
  assert(mcp_interface);
  assert(z);
  assert(F);
  assert(info);
  assert(dparam);
  assert(iparam);

  unsigned n = ((SN_generic_path_env*) mcp_interface->interface_data)->n;
  unsigned nnz = ((SN_generic_path_env*) mcp_interface->interface_data)->nnz;

  Options_Interface *o;
  MCP *m;
  MCP_Termination t;
  Information info_path;
  double *tempZ;
  double *tempF;
  double dnnz;
  int i;

#if defined(USE_OUTPUT_INTERFACE)
  Output_Interface outputInterface =
  {
    NULL,
    messageCB,
    NULL
  };
  Output_SetInterface(&outputInterface);
#else
  /* N.B.: the Output_SetLog call does not work when using a DLL:
   * the IO systems of a .exe and .dll do not automatically interoperate */
  Output_SetLog(stdout);
#endif

  o = Options_Create();
  Path_AddOptions(o);
  Options_Default(o);

  DEBUG_PRINTF("%s: Standalone-C Link\n", Path_Version());

  if (n == 0) {
    fprintf(stdout, "\n ** EXIT - solution found (degenerate model).\n");
    (*info) = MCP_Solved;
    Options_Destroy(o);
    return;
  }

  dnnz = MIN(1.0*nnz, 1.0*n*n);
  if (dnnz > INT_MAX) {
    fprintf(stdout, "\n ** EXIT - model too large.\n");
    (*info) = MCP_Error;
    Options_Destroy(o);
    return;
  }
  nnz = (int) dnnz;

  DEBUG_PRINTF("%d row/cols, %d non-zeros, %3.2f%% dense.\n\n",
          n, nnz, 100.0*nnz/(1.0*n*n));
  nnz++;

  m = MCP_Create(n, nnz);
  MCP_SetInterface(m, mcp_interface);

  Options_Read(o, "path.opt");
  Options_SetDouble(o, "con_tol", SN_get_tolerance(dparam)); /* XXX  */
  DEBUG_OR_VERBOSE(Options_Display(o););

  if (verbose > 0)
  {
    info_path.generate_output = Output_Log | Output_Status | Output_Listing;
  }
  else
  {
    info_path.generate_output = 0;
  }
  info_path.use_start = True;
  info_path.use_basics = True;

  t = Path_Solve(m, &info_path);

  tempZ = MCP_GetX(m);
  tempF = MCP_GetF(m);

  for (i = 0; i < n; ++i) {
    z[i] = tempZ[i];
    F[i] = tempF[i];
  }

  DEBUG_PRINTF("ncp_path :: return code from the solver: %d\n", t);
  if (t == MCP_Solved)
  {
    *info = 0;

    if (verbose > 0)
      printf("PATH : NCP Solved, error\n");
  }
  else if (t == MCP_Error)
  {
    *info = 1;
    if (verbose > 0)
      printf("PATH : Error in the solution.\n");
  }
  else if (t == MCP_NoProgress)
  {
    DEBUG_PRINT("ncp_path :: no progress here\n");
    if (info_path.residual < SN_get_tolerance(dparam))
    {
      *info = 0;
      DEBUG_PRINTF("ncp_path :: no progress, but residual  = %e\n", info_path.residual);
    }
    else
    {
      *info = 1;
      DEBUG_PRINTF("ncp_path :: no progress and residual  = %e\n", info_path.residual);
    }
    *info = 0;
  }
  else
  {
    *info = t;
    if (verbose > 0)
      printf("PATH : Other error: %d\n", t);
  }

  SN_set_residual(dparam, info_path.residual);
  SN_set_nb_iters(iparam, info_path.major_iterations);

  MCP_Destroy(m);
  Options_Destroy(o);
  return;
}

#endif
