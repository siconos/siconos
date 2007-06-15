#include <stdlib.h>
#include <stdio.h>
#include "SiconosApiC.h"

#define NDOF 3
#define NDS  3

#define PA10_NDOF 3

extern void testBouncingBallEDXml();
extern void testBouncingBallTSXml();
extern void testBouncingBallTS();
extern void testBouncingBallED();
extern void testThreeBeadsColumn();

int  main(int argc, char *argv[])
{
  int nbTest;

  if (argc == 1)
  {
    nbTest = 2; /* Default Test */
  }
  else
  {
    nbTest = atoi(argv[1]);

  }


  switch (nbTest)
  {
  case 0:
    testBouncingBallTS();
    break;
  case 1:
    testBouncingBallED();
    break;
  case 2:
    testBouncingBallTSXml();
    break;
  case 3:
    testBouncingBallEDXml();
    break;
  case 4:
    testThreeBeadsColumn();
    break;
  default:
    printf("%s: bad index test \n", argv[0]);
    printf("Usage: %s index \n", argv[0]);
    printf("index values : \n");
    printf("0 - BouncingBallTS\n");
    printf("1 - BouncingBallED\n");
    printf("2 - BouncingBallTSXml\n");
    printf("3 - BouncingBallEDXml\n");
    printf("4 - ThreeBeadsColumn\n");

  }

  exit(0);
}


int error(char *Msg)
{
  printf("sictest error: %s\n", Msg);
  exit(-1);
}


void testBouncingBallEDXml()
{

  FILE *fd;

  int status;
  int k;
  int hasnextevent;
  double plot[4];

  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  status = sicLoadModel("./BallED.xml");

  sicInitSimulation();

  k = 0;

  sicHasNextEvent(&hasnextevent);

  while (hasnextevent == 1)
  {

    sicAdvanceToEvent();
    sicProcessEvents();

    sicHasNextEvent(&hasnextevent);

    plot[0] = k;
    status = sicModelgetQ(&plot[1], 0, 0);


    k++;

    fprintf(fd, "%lf %lf\n", plot[0], plot[1]);

  }

  fclose(fd);

}

void testBouncingBallTSXml()
{

  FILE *fd;

  int status, N;
  int k;
  double plot[4];

  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  status = sicLoadModel("./BallTS.xml");

  sicInitSimulation();

  k = 0;

  sicTimeGetN(&N);

  while (k <= N)
  {

    printf("%d \n", k);

    status = sicSTNextStep();

    status = sicSTComputeOneStep();

    /* store data values */
    plot[0] = k;
    status = sicModelgetQ(&plot[1], 0, 0);

    k++;

    fprintf(fd, "%lf %lf\n", plot[0], plot[1]);
  }

  fclose(fd);

}

void testBouncingBallTS()
{
  int nDof = NDOF;                 /* degrees of freedoom fore beads */
  int idBall, idNSlaw, idRel;      /* id to identify DS, non-smooth law, relation*/
  int  idInter;                    /* id to identify Interaction */
  int  idTime;                     /* id to identify Time dicretisation */
  double t0 = 0;                   /* initial computation time*/
  double T = 10;                   /* final computation time */
  double position_init = 1.0;      /* initial position for lowest bead.*/
  double theta[1]  = {0.5};        /* theta for Moreau integrator */
  double R = 0.1;                  /* Ball radius */
  double m = 1;                    /* Ball mass */
  double g = 9.81;                 /* Gravity */

  /* Dynamical System (Beads) parameters */
  double Mass[NDOF * NDOF] = {m, 0, 0, 0, m, 0, 0, 0, 3. / 5 * m*R * R};

  /*  Dynamical System initial conditions */
  double q0[NDOF] = {position_init, 0, 0};
  double v0[NDOF] = {0, 0, 0};

  double weight[NDOF] = { -m * g, 0, 0};

  /* Interaction parameters */
  double H[NDOF] = {1.0, 0, 0};
  double b[1] = {R};

  /* DS and Interaction set */
  int DS[1];

  /* File data tracing */
  FILE *fd;

  /*   int status,idInter; */
  int k, N, status;
  double plot[4];
  /*   double dH; */

  /* Open File */
  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    error("BouncingBallTS::can not open result.dat");
  }

  /* -- The dynamical system -- */
  idBall = sicLagrangianLinearTIDS(nDof, q0, v0, Mass);
  if (idBall == SIC_ERROR)
    error("BouncingBallTS::sicLagrangianLinearTIDS construction");
  DS[0] = idBall;

  /* -- Set external forces (weight) -- */
  if (sicSetFExt(idBall, weight) == SIC_ERROR)
    error("BouncingBallTS::sicSetFExt construction");

  /* -- Interaction ball-floor -- */
  idNSlaw = sicNewtonImpactNSL(0.9);
  if (idNSlaw == SIC_ERROR)
    error("BouncingBallTS::icNewtonImpactNSL construction");

  H[0] = 1;
  b[0] = R;
  idRel = sicLagrangianLinearR(nDof, 1, H, b);
  if (idRel == SIC_ERROR)
    error("BouncingBallTS::sicLagrangianLinearR construction");


  idInter = sicInteraction("floor-ball", 1, DS, idNSlaw, idRel, 1);
  if (idInter == SIC_ERROR)
    error("BouncingBallTS::sicInteraction construction");

  /* --- NonSmoothDynamicalSystem --- */
  sicNonSmoothDynamicalSystem(0);

  /* --- Model --- */
  sicModel(t0, T);

  /* ------------------ */
  /* --- Simulation --- */
  /* ------------------ */

  idTime = sicTimeDiscretisation(0.001);
  if (idTime == SIC_ERROR)
    error("BouncingBallED::sicTimeDiscretisation construction");

  sicSimulationTimeStepping(idTime);
  sicOneStepIntegratorMoreau(theta);
  sicOneStepNSProblemLCP("Lemke", 101, 0.0001);


  sicInitSimulation();


  k = 0;

  sicTimeGetN(&N);

  while (k < N)
  {

    printf("%d \n", k);

    status = sicSTNextStep();

    status = sicSTComputeOneStep();

    /* store data values */
    plot[0] = k;
    status = sicModelgetQ(&plot[1], 0, 0);

    k++;

    fprintf(fd, "%lf %lf\n", plot[0], plot[1]);
  }

  sicClean();

  fclose(fd);
}

void testBouncingBallED()
{
  int nDof = NDOF;                 /* degrees of freedoom fore beads */
  int idBall, idNSlaw, idRel;      /* id to identify DS, non-smooth law, relation*/
  int  idInter;                    /* id to identify Interaction */
  int  idTime;                     /* id to identify Time dicretisation */
  double t0 = 0;                   /* initial computation time*/
  double T = 10;                   /* final computation time */
  double position_init = 1.0;      /* initial position for lowest bead.*/
  double R = 0.1;                  /* Ball radius */
  double m = 1;                    /* Ball mass */
  double g = 9.81;                 /* Gravity */

  /* Dynamical System (Beads) parameters */
  double Mass[NDOF * NDOF] = {m, 0, 0, 0, m, 0, 0, 0, 3. / 5 * m*R * R};

  /*  Dynamical System initial conditions */
  double q0[NDOF] = {position_init, 0, 0};
  double v0[NDOF] = {0, 0, 0};

  double weight[NDOF] = { -m * g, 0, 0};

  /* Interaction parameters */
  double H[NDOF] = {1.0, 0, 0};
  double b[1] = {R};

  /* DS and Interaction set */
  int DS[1];

  /* File data tracing */
  FILE *fd;

  /*   int status,idInter; */
  int k, status;
  double plot[4];
  int hasnextevent;

  fd = fopen("result.dat", "w");

  if (fd == NULL)
  {
    printf("error:: result.dat write\n");
  }

  /* -- The dynamical system -- */
  idBall = sicLagrangianLinearTIDS(nDof, q0, v0, Mass);
  if (idBall == SIC_ERROR)
    error("BouncingBallTS::sicLagrangianLinearTIDS construction");
  DS[0] = idBall;

  /* -- Set external forces (weight) -- */
  if (sicSetFExt(idBall, weight) == SIC_ERROR)
    error("BouncingBallTS::sicSetFExt construction");

  /* -- Interaction ball-floor -- */
  idNSlaw = sicNewtonImpactNSL(0.9);
  if (idNSlaw == SIC_ERROR)
    error("BouncingBallTS::icNewtonImpactNSL construction");

  H[0] = 1;
  b[0] = R;
  idRel = sicLagrangianLinearR(nDof, 1, H, b);
  if (idRel == SIC_ERROR)
    error("BouncingBallTS::sicLagrangianLinearR construction");


  idInter = sicInteraction("floor-ball", 1, DS, idNSlaw, idRel, 1);
  if (idInter == SIC_ERROR)
    error("BouncingBallTS::sicInteraction construction");

  /* --- NonSmoothDynamicalSystem --- */
  sicNonSmoothDynamicalSystem(0);

  /* --- Model --- */
  sicModel(t0, T);

  /* ------------------ */
  /* --- Simulation --- */
  /* ------------------ */

  idTime = sicTimeDiscretisation(0.001);
  if (idTime == SIC_ERROR)
    error("BouncingBallED::sicTimeDiscretisation construction");

  sicSimulationTimeStepping(idTime);
  sicOneStepIntegratorLsodar(theta);
  sicOneStepNSProblemLCP("LCP", 101, 0.0001);


  sicInitSimulation();

  k = 0;

  sicHasNextEvent(&hasnextevent);

  while (hasnextevent == 1)
  {

    sicAdvanceToEvent();
    sicProcessEvents();

    sicHasNextEvent(&hasnextevent);

    plot[0] = k;
    status = sicModelgetQ(&plot[1], 0, 0);


    k++;

    fprintf(fd, "%lf %lf\n", plot[0], plot[1]);

  }

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



