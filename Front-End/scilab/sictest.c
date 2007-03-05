#include <stdlib.h>
#include <stdio.h>
#include "SiconosApiC.h"

#define NDOF 3
#define NDS  3

#define PA10_NDOF 3

extern void testHuMAns_pa10();
extern void testMultiBeadsColumn();
extern void testThreeBeadsColumn();

main()
{
  testThreeBeadsColumn();
  /* testMultiBeadsColumn();*/
  /*testHuMAns_pa10();*/
  exit(0);
}


int error(char *Msg)
{
  printf("sictest error: %s\n", Msg);
  exit(-1);
}

void testHuMAns_pa10()
{
  /*  Dynamical System initial conditions */
  double q0[PA10_NDOF] = {1, 0, 0};
  double v0[PA10_NDOF] = {0, 0, 0};

  int nDof = 3;
  /* begin and final computation time */
  double t0 = 0, T = 10.0;
  double h = 0.002; /* time step */
  double criterion = 0.001;
  int maxIter = 50;
  int e = 0.9; /* nslaw */
  int e2 = 0.0;

  /* File data tracing */
  FILE *fd;

  /* Simulation variables */

  int status, nId, idInter;
  int k, N;
  double plot[4];
  double dH;
  /* Interaction parameters */
  int DS[1];
  double H[6][3];
  double b[6] = {1.7, 1.7, 0.3, 0.3, 3.14, 3.14};
  double Theta[1];

  /* Dynamical PA10 system creation */
  nId = sicLagrangianDS(nDof, q0, v0);
  /* external plug-in */
  sicSetComputeMassFunction(nId, "RobotPlugin.so", "mass");
  sicSetComputeNNLFunction(nId, "RobotPlugin.so", "NNL");
  sicSetComputeJacobianQNNLFunction(nId, "RobotPlugin.so", "jacobianQNNL");
  sicSetComputeJacobianVelocityNNLFunction(nId, "RobotPlugin.so", "jacobianVNNL");

  sicSetComputeFIntFunction(nId, "RobotPlugin.so", "FInt");
  sicSetComputeJacobianQFIntFunction(nId, "RobotPlugin.so", "jacobianQFInt");
  sicSetComputeJacobianVelocityFIntFunction(nId, "RobotPlugin.so", "jacobianQFInt");
  sicSetComputeFExtFunction(nId, "RobotPlugin.so", "FExt");

  /* -------------------
   * --- Interactions---
   * -------------------

   *  Two interactions:
   *  - one with Lagrangian non linear relation to define contact with ground
   *  - the other to define angles limitations (articular stops), with lagrangian linear relation
   *  Both with newton impact nslaw.
   */

  DS[0] = 0.0;
  idInter = sicInteraction("floor-arm", 1, DS, 2);
  //sicLagrangianR(idInter,"scleronomic", "RobotPlugin:h2","RobotPlugin:G2");
  sicNewtonImpactNSL(idInter, e);

  H[0][0] = -1;
  H[1][0] = 1;
  H[2][1] = -1;
  H[3][1] = 1;
  H[4][2] = -1;
  H[5][2] = 1;

  b[0] = 1.7;
  b[1] = 1.7;
  b[2] = 0.3;
  b[3] = 0.3;
  b[4] = 3.14;
  b[5] = 3.14;
  idInter = sicInteraction("floor-arm2", 1, DS, 6);
  sicLagrangianLinearR(idInter, H, b);
  sicNewtonImpactNSL(idInter, e2);

  /* Construct NSDS */
  sicNonSmoothDynamicalSystem(0);
  /* Construct Model */
  sicModel(t0, T);

  /* Simulation Model */
  sicSimulationTimeStepping(h);
  Theta[0] = 0.5;
  sicOneStepIntegratorMoreau(Theta);
  sicOneStepNSProblemLCP("NLGS", 101, 0.001);

  /* Open File */
  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  /* Simulation */
  sicInitSimulation();

  k = 0;
  sicTimeGetN(&N);

  while (k <= N)
  {
    /* transfer of state i+1 into state i and time incrementation*/
    status = sicSTNextStep();

    /* solve ..*/
    sicSTnewtonSolve(criterion, maxIter);
    /* update */
    status = sicSTupdateState();

    /* --- Get values to be plotted ---*/
    status = sicTimeGetH(&dH);

    plot[0] = k * dH;
    status = sicModelgetQ(&plot[1], 0, 0);
    status = sicModelgetQ(&plot[2], 0, 1);
    status = sicModelgetQ(&plot[3], 0, 2);

    k++;

    fprintf(fd, "%lf %lf %lf %lf\n", plot[0], plot[1], plot[2], plot[3]);
  }

  sicClean();

  fclose(fd);
}


void testMultiBeadsColumn()
{
  int dsNumber = NDS;   /* The number of dynamical systems */
  int nDof = NDOF;      /* degrees of freedoom fore beads */
  double inc_pos = 0.5; /* increment position from one DS to following */
  double inc_vel = 0;   /* increment velocity from one DS to following */
  char nameInter[32];
  int i, index;

  /* Dynamical System (Beads) parameters */
  double Mass[NDOF * NDOF] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double C[NDOF * NDOF] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  double K[NDOF * NDOF] = {0, 0, 0, 0, 1, 0, 0, 0, 0};

  /*  Dynamical System initial conditions */
  double q0[NDOF] = {1, 0, 0};
  double v0[NDOF] = {0, 0, 0};

  /* Interaction parameters */
  int DS[2];
  double H[NDOF * 2] = {0, 0, 0, 0, 0, 0};
  double b[1] = {0};

  /* Simulation parameters */
  double Theta[NDS];

  /* File data tracing */
  FILE *fd;

  /* Simulation variables */

  int status, idInter;
  int k, N;
  double plot[4];
  double dH;

  /* Open File */
  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  /* Creation of Dynamical Systems (Beads) */
  for (i = 0; i < dsNumber; i++)
  {
    if (sicLagrangianLinearTIDS(nDof, q0, v0, Mass, K, C, "BeadsPlugin.so", "beadsFExt") == SIC_ERROR)
      error("sicLagrangianLinearTIDS creation");
    q0[0] += inc_pos;
    v0[0] += inc_vel;
    Theta[i] = 0.500001;
  }

  /* Creation of Interactions */
  /* Bead and floor */
  DS[0] = 0;
  DS[1] = 0;
  H[0] = 1;
  b[0] = 0;
  idInter = sicInteraction("floor", 1, DS, 1);
  sicLagrangianLinearR(idInter, H, b);
  sicNewtonImpactNSL(idInter, 0.9);
  /* Last Bead and ceiling */
  DS[0] = dsNumber - 1;
  DS[1] = 0;
  H[0] = -1;
  b[0] = 1.5;
  idInter = sicInteraction("ceiling", 1, DS, 1);
  sicLagrangianLinearR(idInter, H, b);
  sicNewtonImpactNSL(idInter, 0.9);
  /* Between beads */
  H[0] = -1;
  H[3] = 1;
  b[0] = 0;
  for (i = 1; i < dsNumber; i++)
  {
    DS[0] = i - 1;
    DS[1] = i;
    sprintf(nameInter, "inter%d\0", i);
    idInter = sicInteraction(nameInter, 2, DS, 1);
    sicLagrangianLinearR(idInter, H, b);
    sicNewtonImpactNSL(idInter, 0.9);
  }

  /*Construct NSDS */
  sicNonSmoothDynamicalSystem(0);
  /*Construct Model */
  sicModel(0.0, 10.0);

  /* Simulation Model */
  sicSimulationTimeStepping(0.001);
  sicOneStepIntegratorMoreau(Theta);
  sicOneStepNSProblemLCP("NSQP", 101, 0.0001);


  /* Simulation */
  sicInitSimulation();

  k = 0;
  sicTimeGetN(&N);

  while (k <= N)
  {


    status = sicSTNextStep();

    status = sicSTComputeOneStep();

    status = sicTimeGetH(&dH);

    plot[0] = k * dH;
    index = 0;
    status = sicModelgetQ(&plot[1], index, 0);
    index = 1;
    status = sicModelgetQ(&plot[2], index, 0);
    index = 2;
    status = sicModelgetQ(&plot[3], index, 0);
    k++;

    fprintf(fd, "%lf %lf %lf %lf \n", plot[0], plot[1], plot[2], plot[3]);
  }

  sicClean();

  fclose(fd);
}

void testThreeBeadsColumn()
{
  FILE *fd;

  int status;
  int k, N;
  int index;
  double plot[4];
  double H;

  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  status = sicLoadModel("./ThreeBeadsColumn.xml");

  sicInitSimulation();

  k = 0;
  sicTimeGetN(&N);

  while (k <= N)
  {


    printf("%d \n", k);

    status = sicSTNextStep();

    status = sicSTComputeOneStep();

    status = sicTimeGetH(&H);

    plot[0] = k * H;
    index = 0;
    status = sicModelgetQ(&plot[1], index, 0);
    index = 1;
    status = sicModelgetQ(&plot[2], index, 0);
    index = 2;
    status = sicModelgetQ(&plot[3], index, 0);

    k++;

    fprintf(fd, "%lf %lf %lf %lf \n", plot[0], plot[1], plot[2], plot[3]);

    //sicDebug(&status);

  }

  fclose(fd);
}



