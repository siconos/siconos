/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file FrictionalContact3D.cpp
  A very simple example to chow how to use SiconosNumerics to solve the problem
 * Find \f$(reaction,velocity)\f$ such that:\n
 * \f$
 \left\lbrace
  \begin{array}{l}
  M globalVelocity =  q +  H reaction \\
  velocity = H^T globalVelocity + b\\
  K \ni reaction \perp velocity + \mu \| velocity_t\| \in K^* \\
  \end{array}
  \right.
  \f$\n
  * where
  \f$
  \left\lbrace
  \begin{array}{l}
  K = \{reaction, \|reaction_t\| \leq \mu reaction_n \}
  \end{array}
  \right.
  \f$
  is the Coulomb's Cone \n
    * and with:
    *    - \f$globalVelocity \in R^{n} \f$  the global unknown,
    *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
    *    - \f$velocity \in R^{m} \f$  and \f$reaction \in R^{m} \f$ the local unknowns,
    *    - \f$b \in R^{m} \f$ is the modified local velocity (\f$ e U_{N,k}\f$)
    *    - \f$M \in R^{n \times n } \f$  and \f`$q \in R^{n} \f$
    *    - \f$H \in R^{n \times m } \f$
    \f$ reaction_n\f$ represents the normal part of the reaction while \f$ reaction_t\f$ is its tangential part.

    \f$ \mu \f$ is the friction coefficient (it may be different for each contact).




  \section pfc3DSolversList Available solvers for Friction Contact 3D
  Use the generic function globalFrictionContact3D_driver() to call one the the specific solvers listed below:

  - globalfrictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

  (see the functions/solvers list in GlobalFrictionContact3D_Solvers.h)

  \section pfc3DParam Required and optional parameters
  GlobalFrictionContact3D problems needs some specific parameters, given to the GlobalFrictionContact3D_driver() function thanks to a SolverOptions structure. \n

  \brief
*/
#include "GlobalFrictionContactProblem.h"
#include "FrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "SiconosNumerics.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SolverOptions.h"
#include "Friction_cst.h"
#include "gfc3d_Solvers.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int info = -1;



  int NC = 1;//Number of contacts
  int Ndof = 9;//Number of DOF
  // Problem Definition
  double M11[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M22[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M33[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  /*     double M[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 1, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 1, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 1, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 1, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 1, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 1, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 1, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 0, 1}; */


  double H00[9] =  {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double H20[9] =  { -1, 0, 0, 0, -1, 0, 0, 0, -1};

  /*     double H[27] = {1, 0, 0, 0, 0, 0, -1, 0, 0, */
  /*        0, 1, 0, 0, 0, 0, 0, -1, 0, */
  /*        0, 0, 1, 0, 0, 0, 0, 0, -1}; */


  double q[9] = { -3, -3, -3, -1, 1, 3, -1, 1, 3};
  double b[3] = {0, 0, 0};
  double mu[1] = {0.1};

  /*    DSCAL(9,-1.0,q,1); */




  /*     int NC = 3;//Number of contacts  */
  /*     int Ndof = 9;//Number of DOF  */
  /*     double M[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 1, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 1, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 1, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 1, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 1, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 1, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 1, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 0, 1}; */
  /*     double H[81] = {1, 0, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 1, 0, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 1, 0, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 1, 0, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 1, 0, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 1, 0, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 1, 0, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 1, 0,  */
  /*        0, 0, 0, 0, 0, 0, 0, 0, 1}; */


  /*     double q[9] = {-1, 1, 3, -1, 1, 3, -1, 1, 3}; */
  /*     double b[9] = {0, 0, 0,0, 0, 0,0, 0, 0 }; */
  /*     double mu[3] = {0.1,0.1,0.1};    */


  int k;
  int m = 3 * NC;
  int n = Ndof;

  GlobalFrictionContactProblem numericsProblem;
  globalFrictionContact_null(&numericsProblem);
  numericsProblem.numberOfContacts = NC;
  numericsProblem.dimension = 3;
  numericsProblem.mu = mu;
  numericsProblem.q = q;
  numericsProblem.b = b;

  numericsProblem.M = NM_new();
  NumericsMatrix *MM =  numericsProblem.M;
  MM->storageType = NM_SPARSE_BLOCK;
  MM->size0 = Ndof;
  MM->size1 = Ndof;


  MM->matrix1 = SBM_new();
  SparseBlockStructuredMatrix *MBlockMatrix = MM->matrix1;
  MBlockMatrix->nbblocks = 3;
  double * block[3] = {M11, M22, M33};
  MBlockMatrix->block = block;
  MBlockMatrix->blocknumber0 = 3;
  MBlockMatrix->blocknumber1 = 3;
  unsigned int blocksize[3] = {3, 6, 9} ;
  MBlockMatrix->blocksize0 = blocksize;
  MBlockMatrix->blocksize1 = blocksize;
  MBlockMatrix->filled1 = 4;
  MBlockMatrix->filled2 = 3;
  size_t index1_data[4] = {0, 1, 2, 3} ;
  size_t index2_data[3] = {0, 1, 2} ;
  MBlockMatrix->index1_data =  index1_data;
  MBlockMatrix->index2_data =  index2_data;


  numericsProblem.H = NM_new();
  NumericsMatrix *HH =  numericsProblem.H;
  HH->storageType = 1;
  HH->size0 = Ndof;
  HH->size1 = 3 * NC;

  HH->matrix1 = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));
  SparseBlockStructuredMatrix *HBlockMatrix = HH->matrix1;
  HBlockMatrix->nbblocks = 2;
  double * hblock[3] = {H00, H20};
  HBlockMatrix->block = hblock;
  HBlockMatrix->blocknumber0 = 3;
  HBlockMatrix->blocknumber1 = 1;
  unsigned int blocksize0[3] = {3, 6, 9} ;
  unsigned int blocksize1[1] = {3} ;
  HBlockMatrix->blocksize0 = blocksize0;
  HBlockMatrix->blocksize1 = blocksize1;
  HBlockMatrix->filled1 = 4;
  HBlockMatrix->filled2 = 2;
  size_t hindex1_data[4] = {0, 1, 1, 2} ;
  size_t hindex2_data[3] = {0, 0} ;
  HBlockMatrix->index1_data =  hindex1_data;
  HBlockMatrix->index2_data =  hindex2_data;

  FILE * foutput = fopen("Example_GlobalFrictionContact_SBM.dat", "w");
  globalFrictionContact_printInFile(&numericsProblem,  foutput);
  fclose(foutput);


  // Unknown Declaration

  double *reaction = (double*)calloc(m, sizeof(double));
  double *velocity = (double*)calloc(m, sizeof(double));
  double *globalVelocity = (double*)calloc(n, sizeof(double));
  // Solver Options
  SolverOptions * numerics_solver_options = (SolverOptions *)malloc(sizeof(SolverOptions));
  //    char solvername[10]= "NSGS";

  /*\warning Must be adpated  for future globalFrictionContact3D_setDefaultSolverOptions*/
  gfc3d_setDefaultSolverOptions(numerics_solver_options, SICONOS_GLOBAL_FRICTION_3D_NSGS);
  numerics_solver_options->dparam[0] = 1e-14;
  numerics_solver_options->iparam[0] = 100000;

  //Driver call
  info = gfc3d_driver(&numericsProblem,
		      reaction , velocity, globalVelocity,
		      numerics_solver_options);
  solver_options_delete(numerics_solver_options);

  free(numerics_solver_options);
  // Solver output
  printf("\n");
  for (k = 0 ; k < m; k++) printf("velocity[%i] = %12.8e \t \t reaction[%i] = %12.8e \n ", k, velocity[k], k , reaction[k]);
  for (k = 0 ; k < n; k++) printf("globalVelocity[%i] = %12.8e \t \n ", k, globalVelocity[k]);
  printf("\n");


  free(reaction);
  free(velocity);
  free(globalVelocity);

  //     SBM_free(MM->matrix1);
  //     SBM_free(HH->matrix1);
  free(MM->matrix1);
  MM->matrix1 = NULL;
  free(HH->matrix1);
  HH->matrix1 = NULL;
  NM_free(MM);
  NM_free(HH);
  free(MM);
  free(HH);
  gfc3d_free_workspace(&numericsProblem);


  /*     while (1) sleep(60); */


  return info;


}
