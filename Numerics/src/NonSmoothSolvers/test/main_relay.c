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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "NonSmoothDrivers.h"
#include "LA.h"
#include "blaslapack.h" /* for dlamch */

//#define PATH_DRIVER

#ifdef PATH_DRIVER
const unsigned short int *__ctype_b;
const __int32_t *__ctype_tolower ;
#endif /*PATH_DRIVER*/
/*
******************************************************************************
*/

int test_relay_series(Relay_Problem * problem, int* solversList)
{
  /*


     solversList[i] = 1 => the corresponding solver will be applied to the input problem
     solversList[i] = 0 => solver ignored.

     0: NLGS
     1: Latin


     MIND TO CHANGE totalNBSolver if you add a new solver in the list

  */

  int totalNBSolver = 3;

  int nbSolvers = 0; /* Real number of solvers called in the tests */
  int i, j;
  for (i = 0; i < totalNBSolver; ++i)
  {
    if (solversList[i] == 1)
      nbSolvers++;
  }

  int nonsymmetric;
  int incx = 1, incy = 1;
  double comp, diff;

  int n = problem->size;
  /*   printf ("problem size = %i",n); */
  double **z = malloc(nbSolvers * sizeof(*z));
  double **w = malloc(nbSolvers * sizeof(*w));

  for (i = 0; i < nbSolvers; i++)
  {
    z[i] = malloc(n * sizeof(double));
    w[i] = malloc(n * sizeof(double));

    for (j = 0; j < n; ++j)
    {
      z[i][j] = 0.0;
      w[i][j] = 0.0;
    }
  }
  char nameList[300] = "";

  /* Buffer for w, to check if w = Mz+q  */
  double *wBuffer = malloc(n * sizeof(double));

  Numerics_Options global_options;
  global_options.verboseMode = 0;

  nonsymmetric = 0;
  int info = 0;
  double alpha = -1, beta  = 1;
  int maxIter = 100001;
  double tolerance = 1e-8;

  Solver_Options * options ;
  Solver_Options * local_options = NULL;
  int numberOfSolvers ;
  int isSparse = 0;
  /* Dense storage */
  if (problem->M->storageType == 0)
  {
    printf("\n\n  ");
    printf("The matrix of the RELAY is dense (ie double* storage) ");
    printf("\n\n  ");
    numberOfSolvers = 1;
    options = malloc(numberOfSolvers * sizeof(*options));
    local_options = options;
    /* Is M symmetric ? */
    for (i = 0 ; i < n ; ++i)
    {
      for (j = 0 ; j < i ; ++j)
      {
        if (abs(problem->M->matrix0[i * n + j] - problem->M->matrix0[j * n + i]) > DBL_EPSILON)
        {
          nonsymmetric = 1;
          break;
        }
      }
    }

    if (nonsymmetric) printf("\n !! WARNING !!\n M is a non symmetric matrix \n");
    else printf(" M is a symmetric matrix \n");

  }
  /* Sparse Block storage */
  else
  {
    printf("\n\n  ");
    printf("The matrix of the RELAY is a SparseBlockStructuredMatrix.");
    printf("\n\n  ");
    isSparse = 1;
    numberOfSolvers = 2;
    options = malloc(numberOfSolvers * sizeof(*options));
    strcpy(options[0].solverName, "GaussSeidel_SBM");
    int iparam[3] = {maxIter, 0, 0};
    double dparam[3] = {tolerance, 0.0, 0.0};
    options[0].iSize = 3;
    options[0].dSize = 3;
    options[0].iparam = iparam;
    options[0].dparam = dparam;
    options[0].isSet = 1;
    options[0].filterOn = 0;
    local_options = &options[1];
  }

  printf("\n\n  ");
  printf("The following solvers are called:");
  printf("\n\n  ");
  printf("      SOLVER     | ITER/PIVOT |   ERROR  |  ||w-Mz-q|| |");
  printf("\n\n  ");

  /* Current solver number */
  int k = 0;

  /* NLGS */
  if (solversList[0] == 1)
  {
    strcat(nameList, "    PGS     |");
    strcpy(local_options->solverName, "PGS");
    int iparam[2] = {maxIter, 0};
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 2;
    local_options->dSize = 2;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    local_options->filterOn = 0;
    int info1 = relay_driver(problem, z[k] , w[k], options, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prodNumericsMatrix(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer  , incx);

    printf("  PGS     (LOG:%1d)|      %5d | %10.4g  | %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], diff);

    /*       if(info1!=0) */
    /*  info = info1; */
    k++;
  }

#ifdef HAVE_PATHFERRIS
  if (solversList[2] == 1)
  {
    strcat(nameList, "    PATH     |");
    strcpy(local_options->solverName, "Path");
    double dparam[2] = {tolerance, 0.0};
    local_options->iSize = 0;
    local_options->dSize = 2;
    local_options->iparam = NULL;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    int info1 = relay_driver(problem, z[k] , w[k], options, &global_options);
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer , incy);
    prodNumericsMatrix(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    if (isSparse == 0)
      printf("    Path    (LOG:%1d)|  not set   | %10.4g | %10.4g | %10.4g |\n", info1, local_options->dparam[1], comp, diff);
    else
      printf("    Gauss-Seidel/Path    (LOG:%1d)|  not set   | %10.4g | %10.4g | %10.4g |\n", info1, options[0].dparam[1], comp, diff);

    if (info1 != 0)
      info = info1;
    k++;
  }
#endif
  /* RPGS */
  if (solversList[1] == 1)
  {
    strcat(nameList, "    Latin     |");
    strcpy(local_options->solverName, "Latin");
    int iparam[2] = {maxIter, 0};
    double dparam[3] = {tolerance, 0.0, 1.0};
    local_options->iSize = 2;
    local_options->dSize = 3;
    local_options->iparam = iparam;
    local_options->dparam = dparam;
    local_options->isSet = 1;
    // int info1 = relay_driver( problem, z[k] , w[k], options,  &global_options );
    int info1 = 0;
    comp = DDOT(n , z[k] , incx , w[k] , incy);
    DCOPY(n , w[k], incx, wBuffer , incy);
    DAXPY(n , alpha , problem->q , incx , wBuffer  , incy);
    prodNumericsMatrix(n, n, beta, problem->M, z[k], alpha, wBuffer);
    diff = DNRM2(n , wBuffer , incx);

    printf("    Latin   (LOG:%1d)|      %5d | %10.4g |  %10.4g |\n", info1, local_options->iparam[1], local_options->dparam[1], diff);

    if (info1 != 0)
      /*  info = info1; */
      k++;
  }
  /* =========================== Ouput: comparison between the different methods =========================== */

  strcat(nameList, "\n");
  printf(" *****   z = \n\n");
  printf(nameList);
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < nbSolvers; j++)
      printf("%11.5g | ", z[j][i]);
    printf("\n");
  }
  printf("\n\n");
  printf(" *****   w = \n\n");
  printf(nameList);
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < nbSolvers; j++)
      printf("%11.5g | ", w[j][i]);
    printf("\n");
  }
  printf("\n\n");

  for (i = 0; i < nbSolvers; i++)
  {
    free(z[i]);
    free(w[i]);
  }
  free(z);
  free(w);
  free(wBuffer);
  free(options);
  return info;




}

int test_mmc(void)
{
  printf("========================================================================================================== \n");
  printf("                                     RELAY Solvers tests (function: test_mmc)  \n");
  printf("==========================================================================================================\n");
  FILE *f1, *f2, *f3, *f4;
  int i, nl, nc;
  double qi, Mij;
  char val[20];

  /* Building of the RELAY */
  Relay_Problem * problem = malloc(sizeof(*problem));

  // computation of the size of the problem

  int n = 0;
  // Open files ...

  if ((f1 = fopen("DATA/M_relay1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }

  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    n = nl;
  }

  printf("\n SIZE OF THE PROBLEM : %d \n", n);
  fclose(f1);


  problem->size = n;


  // M
  double * vecM = (double*)malloc(n * n * sizeof(double));
  double * q       = (double*) malloc(n * sizeof(double));
  double * a       = (double*) malloc(n * sizeof(double));
  double * b       = (double*) malloc(n * sizeof(double));

  for (i = 0 ; i < n * n ; ++i)
  {
    vecM[i] = 0.0;
  }
  for (i = 0 ; i < n ; ++i)
  {
    q[i] = 0;
    a[i] = 0;
    b[i] = 0;
  }

  /* Data loading for M and q */

  if ((f1 = fopen("DATA/M_relay1.dat", "r")) == NULL)
  {
    perror("fopen 1");
    exit(1);
  }


  if ((f2 = fopen("DATA/q_relay1.dat", "r")) == NULL)
  {
    perror("fopen 2");
    exit(2);
  }


  if ((f3 = fopen("DATA/a_relay1.dat", "r")) == NULL)
  {
    perror("fopen 5");
    exit(5);
  }


  if ((f4 = fopen("DATA/b_relay1.dat", "r")) == NULL)
  {
    perror("fopen 6");
    exit(6);
  }



  while (!feof(f1))
  {
    fscanf(f1, "%d", &nl);
    fscanf(f1, "%d", &nc);
    fscanf(f1, "%s", val);
    Mij = atof(val);


    vecM[(nc - 1)*n + nl - 1 ] = Mij;

  }




  while (!feof(f2))
  {
    fscanf(f2, "%d", &nl);
    fscanf(f2, "%s", val);
    qi = atof(val);
    q[nl - 1] = -qi;

  }


  while (!feof(f3))
  {
    fscanf(f3, "%d", &nl);
    fscanf(f3, "%s", val);
    qi = atof(val);
    a[nl - 1] = qi;
  }


  while (!feof(f4))
  {
    fscanf(f4, "%d", &nl);
    fscanf(f4, "%s", val);
    qi = atof(val);
    b[nl - 1] = qi;
  }


  fclose(f2);
  fclose(f4);
  fclose(f1);
  fclose(f3);



  NumericsMatrix * MM = malloc(sizeof(*MM));
  MM->matrix0 = vecM;
  MM->size0 = n;
  MM->size1 = n;
  MM->storageType = 0;



  // Set M and q of the problem
  problem->q = q;
  problem->M = MM;
  problem->ub = a;
  problem->lb = b;

  // Call tests

  printf(" ----------------------------------------------------------\n");
  printf("Run working tests ...\n");
  /* Stable: */
  //  int solversList[11] ={1,1,1,1,0,0,0,1,1,0,0};
  int solversList[3] = {1, 1, 1};
  int info = test_relay_series(problem, solversList);
  printf(" ----------------------------------------------------------\n");

  /*   /\* Fail or unstable: *\/ */
  /*   printf("---------------------------------------------------------- \n"); */
  /*   printf("\n Run unstable tests (results may be wrong or log !=0)...\n"); */
  /*   int solversList2[12] ={0,0,0,0,0,0,0,0,0,0,1,1}; */
  /*   int infoFail = test_relay_series(problem, solversList2); */
  /*   printf("--------- End of unstable tests --------------------------- \n"); */

  // Release memory
  problem->M = NULL;
  problem->q = NULL;
  problem->ub = NULL;
  problem->lb = NULL;
  free(MM);
  free(problem);
  free(vecM);
  free(q);
  free(a);
  free(b);

  printf("========================================================================================================== \n");
  printf("                                           END OF TEST MMC     \n");
  printf("==========================================================================================================\n");
  return info;

}

/* To read in a file a Relay_Problem with a "double*" storage for M */
void getProblem(char* name, Relay_Problem *  problem)
{

  FILE * Relayfile =  fopen(name, "r");
  if (Relayfile == NULL)
  {
    fprintf(stderr, "fopen Relayfile: %s\n", name);
    exit(1);
  }
  printf("\n\n******************************************************\n");
  printf("Read Relay Complementarity Problem in file %s\n", name);
  printf("******************************************************\n");

  /* Dim of the Relay */
  int dim;
  fscanf(Relayfile , "%d" , &dim);
  int dim2 = dim * dim;

  problem->M->matrix0 = malloc(dim2 * sizeof(double));
  problem->q = (double*)malloc(dim * sizeof(double));
  problem->lb = (double*) malloc(dim * sizeof(double));
  problem->ub = (double*) malloc(dim * sizeof(double));
  double * vecM = problem->M->matrix0;

  int i, j;
  char val[20];
  /* fill M */
  for (i = 0 ; i < dim ; ++i)
  {
    for (j = 0 ; j < dim ; ++j)
    {
      fscanf(Relayfile, "%s", val);
      vecM[ dim * j + i ] = atof(val);
    }
  }

  /* fill q */
  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(Relayfile , "%s" , val);
    problem->q[i] = atof(val);
  }
  /* fill lb */
  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(Relayfile , "%s" , val);
    problem->lb[i] = atof(val);
  }
  /* fill ub */
  for (i = 0 ; i < dim ; ++i)
  {
    fscanf(Relayfile , "%s" , val);
    problem->ub[i] = atof(val);
  }

  /* fill sol */
  double* sol = NULL;
  fscanf(Relayfile , "%s" , val);
  if (!feof(Relayfile))
  {
    sol  = (double*)malloc(dim * sizeof(double));
    sol[0] = atof(val);
    for (i = 1 ; i < dim ; ++i)
    {
      fscanf(Relayfile , "%s" , val);
      sol[i] = atof(val);
    }
  }
  printf("\n exact solution : ");
  if (sol != NULL) for (i = 0 ; i < problem->size ; ++i) printf(" %10.4g " , sol[i]);
  else printf(" unknown ");
  printf("\n");

  problem->size = dim;
  problem->M->size0 = dim;
  problem->M->size1 = dim;

  fclose(Relayfile);
  if (sol != NULL)
    free(sol);
}

/* To read in a file a LinearComplementarity_Problem with a "SparseBlockStructuredMatrix*" storage for M */
void getProblemSBM(char* name, Relay_Problem *  problem)
{

  /*   FILE * RELAYfile =  fopen( name,"r" ); */
  /*   if(RELAYfile==NULL) */
  /*     { */
  /*       fprintf(stderr,"fopen RELAYfile: %s\n",name); */
  /*       exit(1); */
  /*     } */
  /*   printf("\n\n******************************************************\n"); */
  /*   printf("Read Linear Complementarity Problem in file %s\n", name); */
  /*   printf("******************************************************\n"); */

  /*   printf("\n The matrix M of the RELAY is a SparseBlockStructuredMatrix.\n"); */

  /*   SparseBlockStructuredMatrix * blmat =  problem->M->matrix1; */

  /*   int i,j; */
  /*   char val[20]; */

  /*   /\***** M *****\/ */
  /*   fscanf( RELAYfile , "%d" , &blmat->nbblocks ); */
  /*   fscanf( RELAYfile , "%d" , &blmat->size ); */
  /*   blmat->blocksize = (int*)malloc( blmat->size * sizeof(int) ); */
  /*   for( i = 0 ; i < blmat->size ; i++) fscanf( RELAYfile , "%d" , &blmat->blocksize[i] ); */
  /*   blmat->RowIndex = (int*)malloc( blmat->nbblocks * sizeof(int) ); */
  /*   blmat->ColumnIndex = (int*)malloc( blmat->nbblocks * sizeof(int) ); */
  /*   for( i = 0 ; i < blmat->nbblocks ; i++ ) */
  /*     { */
  /*       fscanf( RELAYfile , "%d" , &blmat->RowIndex[i] ); */
  /*       fscanf( RELAYfile , "%d" , &blmat->ColumnIndex[i] ); */
  /*     } */

  /*   blmat->block = (double**)malloc( blmat->nbblocks * sizeof(double*) ); */
  /*   int pos, sizebl, numberOfRows, numberOfColumns; */
  /*   for (i = 0 ; i < blmat->nbblocks ; i++)  */
  /*     { */
  /*       pos = blmat->RowIndex[i]; */
  /*       numberOfRows = blmat->blocksize[pos]; */
  /*       if(pos>0) */
  /*  numberOfRows -= blmat->blocksize[pos-1]; */
  /*       pos = blmat->ColumnIndex[i]; */
  /*       numberOfColumns = blmat->blocksize[pos]; */
  /*       if(pos>0) */
  /*  numberOfColumns -= blmat->blocksize[pos-1]; */
  /*       sizebl = numberOfRows * numberOfColumns; */
  /*       blmat->block[i] = (double*)malloc( sizebl * sizeof(double) ); */
  /*       for( j = 0 ; j < sizebl ; j++ ){ */
  /*  fscanf(RELAYfile,"%s",val); */
  /*  blmat->block[i][j] = atof(val); */
  /*       } */
  /*     } */

  /*   int dim = blmat->blocksize[blmat->size-1]; */
  /*   /\**** q ****\/ */
  /*   problem->q = (double*)malloc( dim * sizeof(double) ); */
  /*   for( i = 0 ; i < dim ; i++ ){ */
  /*     fscanf( RELAYfile , "%s" , val ); */
  /*     problem->q[i] = atof( val ); */
  /*   } */

  /*   fscanf( RELAYfile , "%s" , val ); */

  /*   double* sol = NULL; */
  /*   if( !feof( RELAYfile ) ){ */
  /*     sol  = (double*)malloc(  dim*sizeof( double ) ); */
  /*     sol[0] = atof( val ); */
  /*     for( i = 1 ; i < dim ; i++ ){ */
  /*       fscanf( RELAYfile , "%s" , val ); */
  /*       sol[i] = atof( val ); */
  /*     } */
  /*   } */
  /*   printf("\n exact solution : "); */
  /*   if(sol!=NULL) for( i = 0 ; i < problem->size ; ++i ) printf(" %10.4g " , sol[i] ); */
  /*   else printf(" unknown "); */
  /*   printf("\n"); */

  /*   problem->size = dim; */
  /*   problem->M->size0 = dim; */
  /*   problem->M->size1 = dim; */
  /*   fclose(RELAYfile); */
  /*   if(sol!=NULL) */
  /*     free(sol); */
}

int test_matrix(void)
{
  /*
    Two problems: one with dense (double*) storage for M, "problem", and the other
    with SparseBlockStructuredMatrix storage, "problemSBM".

    For each case, one or both problems are read in a dat file. \n
    The according to the value of the lists solversList and solversList2 (for problem) ,\n
    solversListSBM (for problemSBM).
    the problems are solved with different solvers. \n
    list[i] = 1 => solver is called, list[i] = 0, solver is ignored.
    Check on top of test_relay_series() function for the list of available solvers and their corresponding indices.

    "stable" tests, that must succeed, are those defined in solversList, while "unstable", that may fail or
    return an unexpected termination value are those defined in solversList2.

   */

  printf("========================================================================================================== \n");
  printf("                         RELAY Solvers tests (function: test_matrix)  \n");
  printf("==========================================================================================================\n");

  FILE *RELAYfile = NULL, *RELAYfileBlock = NULL;
  int i, j, itest;
  int iter = 0;
  double criteria = 0.0;
  double *sol = NULL;

  int NBTEST = 1;

  /* === Building of the RELAYs === */

  /* RELAY with dense storage */
  Relay_Problem * problem = malloc(sizeof(*problem));
  problem->M = malloc(sizeof(*(problem->M)));
  problem->M->storageType = 0;
  problem->M->matrix1 = NULL;

  /*   /\* RELAY with sparse-block storage *\/ */
  /*   Relay_Problem * problemSBM = malloc(sizeof(*problemSBM)); */
  /*   problemSBM->M = malloc(sizeof(*(problemSBM->M))); */
  /*   problemSBM->M->storageType = 1; */
  /*   problemSBM->M->matrix0 = NULL; */
  /*   problemSBM->M->matrix1 = malloc(sizeof(*( problemSBM->M->matrix1))); */

  /* List of working solvers */
  int * solversList = NULL; /* for dense */
  int * solversListSBM = NULL; /* for sparse */

  /* List of unstable solvers (failed or wrong termination value) */
  int * solversList2 = NULL;
  int hasSBM = 0; /* =1 if the read matrix exists also in sparse */
  int hasDense = 0;
  int hasUnstable = 0;
  int info = -1;
  for (itest = 0 ; itest < NBTEST ; ++itest)
  {
    hasSBM = 0;
    hasUnstable = 0;
    hasDense = 0;
    /* Read the RELAY */
    switch (itest)
    {
    case 0:
      getProblem("DATA/relay_deudeu.dat", problem);
      hasSBM = 0;
      hasDense = 1;
      hasUnstable = 0;
      {
        int l1[3] = {1, 1, 1};
        solversList = l1;
      }
      break;
    case 1:
      getProblem("MATRIX/trivial.dat", problem);
      /*  getProblemSBM("MATRIX/trivial_block.dat", problemSBM); */
      hasSBM = 1;
      hasDense = 1;
      hasUnstable = 1;
      {
        int l1[12] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0};
        int l2[12] = {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        int l3[12] = {1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0};
        solversList = l1;
        solversList2 = l2;
        solversListSBM = l3;
      }
      break;

    }

    printf(" ----------------------------------------------------------\n");
    printf("Run working tests ...\n");
    /* Stable: */
    int infoTmp = -1;
    if (hasDense == 1)
      infoTmp = test_relay_series(problem, solversList);

    /*       if(hasSBM == 1) */
    /*  { */
    /*    printf("Run working tests for sparse storage ...\n"); */
    /*    infoTmp = test_relay_series(problemSBM, solversListSBM); */
    /*  } */

    printf(" ----------------------------------------------------------\n\n");

    if (hasUnstable == 1)
    {
      /* Fail or unstable: */
      printf("---------------------------------------------------------- \n");
      printf("\n Run unstable tests (results may be wrong or log !=0)...\n");
      int infoFail = test_relay_series(problem, solversList2);
      printf("--------- End of unstable tests --------------------------- \n\n");
    }

    /* Free Memory */
    if (problem->M->matrix0 != NULL)
      free(problem->M->matrix0);
    problem->M->matrix0 = NULL;
    if (problem->q != NULL)
      free(problem->q);
    problem->q = NULL;
    if (problem->lb != NULL)
      free(problem->lb);
    problem->lb = NULL;
    if (problem->ub != NULL)
      free(problem->ub);
    problem->ub = NULL;
    /*     if(problemSBM->q!=NULL) */
    /*  free(problemSBM->q); */
    /*       problemSBM->q = NULL; */
    info = infoTmp;
    /*      if(hasSBM==1) */
    /*  freeSBM(problemSBM->M->matrix1); */

    if (infoTmp != 0)
      break;
  }
  /*   free(problemSBM->M->matrix1); */
  /*   free(problemSBM); */
  free(problem->M);
  free(problem);
  printf("========================================================================================================== \n");
  printf("                                  END OF TEST MATRIX   \n");
  printf("========================================================================================================== \n");

  return info;

}


int main(void)
{

  /* In each test function, two series of tests are called:
     - working tests => must return info = 0 , else error, test fail.
     - unstable tests => just print results on screen

  */

  int info1  = test_mmc();

  if (info1 != 0)
  {
    // printf("Warning: test_mmc log output different from 0 for some solvers.\n");
    return info1;
  }
  int info2 = test_matrix();
  if (info2 != 0)
  {
    // printf("Warning: test_matrix log output different from 0 for some solvers.\n");
    return info2;
  }
  return 0;
}

