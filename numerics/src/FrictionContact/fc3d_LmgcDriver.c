#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "SiconosBlas.h"
#include "fc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "SparseMatrix_internal.h"

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

#include "fclib_interface.h"
#include "numerics_verbose.h"
#include "SiconosCompat.h"
static int fccounter = -1;

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"


int fc3d_LmgcDriver(double *reaction,
                    double *velocity,
                    double *q,
                    double *mu,
                    double* W,
                    unsigned int *row,
                    unsigned int *column,
                    unsigned int nc,
                    unsigned int nb,
                    int solver_id,
                    double tolerance,
                    int itermax,
                    int verbose,
                    int outputFile,
                    int freq_output,
                    int ndof)
{

  numerics_set_verbose(verbose);

  SparseBlockCoordinateMatrix* MC =  SBCM_new_3x3(nc, nc, nb, row, column, W);

  SparseBlockStructuredMatrix* M = SBCM_to_SBM(MC);

  NumericsMatrix* NM = NM_new_SBM(nc * 3, nc * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new(3, nc, NM, q, mu);

  /* frictionContact_display(FC); */

  SolverOptions numerics_solver_options;

  fc3d_setDefaultSolverOptions(&numerics_solver_options, solver_id);

  if (solver_id == SICONOS_FRICTION_3D_NSGS)
  {
    numerics_solver_options.iparam[SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  }
  else if  (solver_id == SICONOS_FRICTION_3D_NSN_AC)
  {
    numerics_solver_options.iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_NO;
    numerics_solver_options.iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY]=SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN;
  }


  numerics_solver_options.dparam[0] = tolerance;
  numerics_solver_options.iparam[0] = itermax;

  double * reaction_guess;
  double * velocity_guess;
  if (outputFile == 3)
  {
    // Save guesses.

    reaction_guess = (double *)malloc(nc*3*sizeof(double));
    velocity_guess = (double *)malloc(nc*3*sizeof(double));
    for (unsigned int k =0; k < 3*nc; k++) reaction_guess[k]=reaction[k];
    for (unsigned int k =0; k < 3*nc; k++) velocity_guess[k]=velocity[k];

  }

  DEBUG_EXPR(frictionContact_display(FC););

  int info = fc3d_driver(FC, reaction , velocity, &numerics_solver_options);


//  uncomment to save FrictionContactProblem

  if (outputFile == 1)
  {
    FILE * file = fopen("tutu.c", "w");

    fprintf(file, "int nc = %i ;\n ", nc);
    fprintf(file, "int nb = %i ;\n ", nb);
    fprintf(file, "double mu[%i] ={\n", nc);
    for (unsigned int i = 0; i < nc - 1 ; i++)
    {
      fprintf(file, "%32.24e, \t", mu[i]);
    }
    fprintf(file, "%32.24e };\n", mu[nc - 1]);
    fprintf(file, "int row[%i] ={\n", nb);
    for (unsigned int i = 0; i < nb - 1 ; i++)
    {
      fprintf(file, "%i,\t", row[i]);
    }

    fprintf(file, " %i};\n", row[nb - 1]);
    fprintf(file, "int column[%i] ={\n", nb);
    for (unsigned int i = 0; i < nb - 1 ; i++)
    {
      fprintf(file, "%i,\t", column[i]);
    }
    fprintf(file, " %i};\n", column[nb - 1]);
    fprintf(file, "double q[%i] ={\n", 3 * nc);
    for (unsigned int i = 0; i < 3 * nc - 1 ; i++)
    {
      fprintf(file, "%32.24e,\t", q[i]);
    }
    fprintf(file, " %32.24e};\n", q[3 * nc - 1]);

    fprintf(file, "double W[%i] ={\n", 3 * 3 * nb);
    for (unsigned int i = 0; i < nb - 1 ; i++)
    {
      for (unsigned int j = 0; j < 3 * 3 ; j++)
      {
        fprintf(file, "%32.24e, \t", W[i * 9 + j]);
      }
      fprintf(file, "\n");
    }
    for (int j = 0; j < 3 * 3 - 1 ; j++)
    {
      fprintf(file, "%32.24e, \t", W[(nb - 1) * 9 + j]);
    }
    fprintf(file, "%32.24e};\n", W[(nb - 1) * 9 + 8]);
    fclose(file);
  }
  else if (outputFile == 2)
  {
    char fname[256];
    sprintf(fname, "LMGC_FC3D-i%.5d-%i-%.5d.dat", numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE], nc, fccounter++);
    printf("LMGC_FC3D-i%.5d-%i-%.5d.dat", numerics_solver_options.iparam[7], nc, fccounter++);
    FILE * foutput  =  fopen(fname, "w");
    frictionContact_printInFile(FC, foutput);
    fclose(foutput);
  }
  else if (outputFile == 3)
  {
#ifdef WITH_FCLIB
    fccounter ++;
    if (fccounter % freq_output == 0)
    {
      char fname[256];
      sprintf(fname, "LMGC_FC3D-i%.5d-%i-%.5d.hdf5", numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE], nc, fccounter);
      printf("Dump LMGC_FC3D-i%.5d-%i-%.5d.hdf5.\n", numerics_solver_options.iparam[SICONOS_IPARAM_ITER_DONE], nc, fccounter);
      /* printf("ndof = %i.\n", ndof); */

      FILE * foutput  =  fopen(fname, "w");
      int n = 100;



      char * title = (char *)malloc(n * sizeof(char));

      strcpy(title, "LMGC dump in hdf5");




      char * description = (char *)malloc(n * sizeof(char));
      strcpy(description, "Rewriting in hdf5 through siconos of  ");
      strcat(description, fname);
      strcat(description, " in FCLIB format");
      char * mathInfo = (char *)malloc(n * sizeof(char));
      strcpy(mathInfo,  "unknown");

      frictionContact_fclib_write(FC,
                                  title,
                                  description,
                                  mathInfo,
                                  fname,ndof);

      frictionContact_fclib_write_guess(
        reaction_guess,
        velocity_guess,
        fname);

      fclose(foutput);
    }
#else
    printf("Fclib is not available ...\n");
#endif
    free(reaction_guess);
    free(velocity_guess);
  }







   SBCM_free_3x3(MC);

  free(M->index1_data);
  free(M->index2_data);
  free(M->block);
  free(M);
  free(FC);
  solver_options_delete(&numerics_solver_options);
  free(NM);
  free(MC);

  return info;
}
