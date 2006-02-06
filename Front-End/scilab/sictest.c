#include <stdio.h>
#include "SiconosApiC.h"

#define NDOF 3
#define NDS  3

main()
{
  /* testThreeBeadsColumn();*/
  testMultiBeadsColumn();
}

int error(char *Msg)
{
  printf("sictest error: %s\n", Msg);
  exit(-1);
}

testMultiBeadsColumn()
{
  int dsNumber = NDS;   /* The number of dynamical systems */
  int nDof = NDOF;      /* degrees of freedoom fore beads */
  double inc_pos = 0.5; /* increment position from one DS to following */
  double inc_vel = 0;   /* increment velocity from one DS to following */
  char nameInter[32];
  int i;

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

  /* Strategy parameters */
  double Theta[NDS];

  /* File data tracing */
  FILE *fd;

  /* Simulation variables */

  int status, idInter;
  int k, N;
  int index;
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
  sicNewtonImpactLawNSL(idInter, "NewtonImpactLawNSL", 0.9);
  /* Last Bead and ceiling */
  DS[0] = dsNumber - 1;
  DS[1] = 0;
  H[0] = -1;
  b[0] = 1.5;
  idInter = sicInteraction("ceiling", 1, DS, 1);
  sicLagrangianLinearR(idInter, H, b);
  sicNewtonImpactLawNSL(idInter, "NewtonImpactLawNSL", 0.9);
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
    sicNewtonImpactLawNSL(idInter, "NewtonImpactLawNSL", 0.9);
  }

  /*Construct NSDS */
  sicNonSmoothDynamicalSystem(0);
  /*Construct Model */
  sicModel(0.0, 10.0);

  /* Strategy Model */
  sicStrategyTimeStepping(0.001);
  sicOneStepIntegratorMoreau(Theta);
  sicOneStepNSProblemLCP(101, 0.0001);


  /* Simulation */
  sicInitStrategy();

  sicTimeGetK(&k);
  sicTimeGetN(&N);

  while (k <= N)
  {


    status = sicSTNextStep();

    status = sicTimeGetK(&k);

    status = sicSTComputeFreeState();

    status = sicSTComputePb();

    status = sicSTupdateState();

    status = sicTimeGetH(&dH);

    plot[0] = k * dH;
    index = 0;
    status = sicModelgetQ(&plot[1], index, 0);
    index = 1;
    status = sicModelgetQ(&plot[2], index, 0);
    index = 2;
    status = sicModelgetQ(&plot[3], index, 0);

    fprintf(fd, "%lf %lf %lf %lf \n", plot[0], plot[1], plot[2], plot[3]);
  }

  sicClean();

  fclose(fd);
}

testThreeBeadsColumn()
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

  sicInitStrategy();

  sicTimeGetK(&k);
  sicTimeGetN(&N);

  while (k <= N)
  {


    printf("%d \n", k);

    status = sicSTNextStep();

    status = sicTimeGetK(&k);

    status = sicSTComputeFreeState();

    status = sicSTComputePb();

    status = sicSTupdateState();

    status = sicTimeGetH(&H);

    plot[0] = k * H;
    index = 0;
    status = sicModelgetQ(&plot[1], index, 0);
    index = 1;
    status = sicModelgetQ(&plot[2], index, 0);
    index = 2;
    status = sicModelgetQ(&plot[3], index, 0);

    fprintf(fd, "%lf %lf %lf %lf \n", plot[0], plot[1], plot[2], plot[3]);

    //sicDebug(&status);

  }

  fclose(fd);
}



