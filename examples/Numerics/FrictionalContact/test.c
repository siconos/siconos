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
  Use the generic function primalFrictionContact3D_driver() to call one the the specific solvers listed below:

  - primalfrictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

  (see the functions/solvers list in GlobalFrictionContact3D_Solvers.h)

  \section pfc3DParam Required and optional parameters
  GlobalFrictionContact3D problems needs some specific parameters, given to the GlobalFrictionContact3D_driver() function thanks to a SolverOptions structure. \n

  \brief
*/

#include "SiconosNumerics.h"

int main(int argc, char* argv[])
{


  // Problem Definition
  int info = -1;



  int NC = 2;//Number of contacts
  int Ndof = 12;//Number of DOF
  // Problem Definition
  double M11[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M22[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M33[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
  double M44[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Warning Fortran Storage
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
  double H11[9] =  {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double H31[9] =  { -1, 0, 0, 0, -1, 0, 0, 0, -1};

  /*     double H[27] = {1, 0, 0, 0, 0, 0, -1, 0, 0, */
  /*        0, 1, 0, 0, 0, 0, 0, -1, 0, */
  /*        0, 0, 1, 0, 0, 0, 0, 0, -1}; */


  double q[12] = { -3, -3, -3, -3, -3, -3, -1, 1, 3, -10, 1, 3};
  double b[6] = {0, 0, 0, 0, 0, 0};
  double mu[2] = {0.1, 0.1};

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


  int i, j, k;
  int m = 3 * NC;
  int n = Ndof;

  GlobalFrictionContactProblem NumericsProblem;
  NumericsProblem.numberOfContacts = NC;
  NumericsProblem.dimension = 3;
  NumericsProblem.mu = mu;
  NumericsProblem.q = q;
  NumericsProblem.b = b;

  NumericsProblem.M = NM_new();
  NumericsMatrix *MM =  NumericsProblem.M;
  MM->storageType = 1;
  MM->size0 = Ndof;
  MM->size1 = Ndof;


  MM->matrix1 = SBM_new();
  MM->matrix0 = NULL;
  SparseBlockStructuredMatrix *MBlockMatrix = MM->matrix1;
  MBlockMatrix->nbblocks = 4;
  double * block[4] = {M11, M22, M33, M44};
  MBlockMatrix->block = block;
  MBlockMatrix->blocknumber0 = 4;
  MBlockMatrix->blocknumber1 = 4;
  int blocksize[4] = {3, 6, 9, 12} ;
  MBlockMatrix->blocksize0 = blocksize;
  MBlockMatrix->blocksize1 = blocksize;
  MBlockMatrix->filled1 = 5;
  MBlockMatrix->filled2 = 4;
  size_t index1_data[5] = {0, 1, 2, 3, 4} ;
  size_t index2_data[4] = {0, 1, 2, 3} ;
  MBlockMatrix->index1_data =  index1_data;
  MBlockMatrix->index2_data =  index2_data;


  NumericsProblem.H = SBM_new();
  NumericsMatrix *HH =  NumericsProblem.H;
  HH->storageType = 1;
  HH->size0 = Ndof;
  HH->size1 = 3 * NC;

  HH->matrix1 = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));
  HH->matrix0 = NULL;
  SparseBlockStructuredMatrix *HBlockMatrix = HH->matrix1;
  HBlockMatrix->nbblocks = 4;
  double * hblock[4] = {H00, H11, H20, H31};
  HBlockMatrix->block = hblock;
  HBlockMatrix->blocknumber0 = 4;
  HBlockMatrix->blocknumber1 = 2;
  int blocksize0[4] = {3, 6, 9, 12} ;
  int blocksize1[3] = {3, 6, 9} ;
  HBlockMatrix->blocksize0 = blocksize0;
  HBlockMatrix->blocksize1 = blocksize1;
  HBlockMatrix->filled1 = 5;
  HBlockMatrix->filled2 = 4;
  size_t hindex1_data[5] = {0, 1, 2, 3, 4} ;
  size_t hindex2_data[4] = {0, 1, 0, 1} ;
  HBlockMatrix->index1_data =  hindex1_data;
  HBlockMatrix->index2_data =  hindex2_data;
  SBM_print(HBlockMatrix);


  // Unknown Declaration

  double *reaction = (double*)malloc(m * sizeof(double));
  double *velocity = (double*)malloc(m * sizeof(double));
  double *globalVelocity = (double*)malloc(n * sizeof(double));
  for (k = 0 ; k < m; k++)
  {
    velocity[k] = 0.0;
    reaction[k] = 0.0;
  }
  for (k = 0 ; k < n; k++) globalVelocity[k] = 0.0;
  // Numerics and Solver Options

  SolverOzptions numerics_solver_options;
  numerics_solver_options.filterOn = 0;
  numerics_solver_options.isSet = 1;

  //    strcpy(numerics_solver_optionslverName,"NSGS_WR");
  //    strcpy(numerics_solver_optionslverName,"NSGS");
  //    strcpy(numerics_solver_optionslverName,"NSGSV_WR");
  numerics_solver_optionslverId = SICONOS_FRICTION_3D_NSGS;

  numerics_solver_options.iSize = 5;
  numerics_solver_options.iparam = (int*)malloc(numerics_solver_options.iSize * sizeof(int));
  int ii ;
  for (ii = 0; ii < numerics_solver_options.iSize; ii++)
    numerics_solver_options.iparam[ii] = 0;
  numerics_solver_options.dSize = 5;
  numerics_solver_options.dparam = (double*)malloc(numerics_solver_options.dSize * sizeof(double));
  for (ii = 0; ii < numerics_solver_options.dSize; ii++)
    numerics_solver_options.dparam[ii] = 0;

  int nmax = 100; // Max number of iteration
  int localsolver = 0; // 0: projection on Cone, 1: Newton/AlartCurnier,  2: projection on Cone with local iteration, 2: projection on Disk  with diagonalization,
  double tolerance = 1e-10;
  double localtolerance = 1e-12;



  numerics_solver_options.iparam[0] = nmax ;
  numerics_solver_options.iparam[4] = localsolver ;
  numerics_solver_options.dparam[0] = tolerance ;
  numerics_solver_options.dparam[2] = localtolerance ;

  //Driver call
  info = primalFrictionContact3D_driver(&NumericsProblem,
                                        reaction , velocity, globalVelocity,
                                        &numerics_solver_options);



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
  free(HH->matrix1);
  free(MM);
  free(HH);

  free(numerics_solver_options.iparam);
  free(numerics_solver_options.dparam);

  /*     while (1) sleep(60); */


  return info;


}
