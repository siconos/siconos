#include "SiconosApiSci.h"


int sicLoadModelInterface(char *fname)
{
  static int minrhs = 1, maxrhs = 1, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int dimo1 = 1, dimo2 = 1, st;
  static int ModelXmlFile;

  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get ModelXmlFile (1, char *)  */
  GetRhsVar(1, "c", &dim1, &dim2, &ModelXmlFile);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLoadModel (number ModelXmlFile has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicLoadModel(cstk(ModelXmlFile));
  /*  Return variables  */
  LhsVar(1) = 2;

  return 0;

}

int sicInitStrategyInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicInitStrategy();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicTimeGetHInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 2;
  static int dim1 = 1, dim2 = 1, H;
  static int dimo1 = 1, dimo2 = 1, st;


  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "d", &dim1, &dim2, &H);
  CreateVar(2, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicTimeGetH(stk(H));
  /* Return variable*/
  LhsVar(1) = 1;
  LhsVar(2) = 2;

  return 0;
}

int sicTimeGetNInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 2;
  static int dim1 = 1, dim2 = 1, N;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dim1, &dim2, &N);
  CreateVar(2, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicTimeGetN(istk(N));
  /* Return variable*/
  LhsVar(1) = 1;
  LhsVar(2) = 2;
  return 0;
}

int sicTimeGetKInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 2;
  static int dim1 = 1, dim2 = 1, K;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dim1, &dim2, &K);
  CreateVar(2, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicTimeGetK(istk(K));

  /* Return variable*/
  LhsVar(1) = 1;
  LhsVar(2) = 2;

  return 0;
}

int sicSTNextStepInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTNextStep();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicSTComputeFreeStateInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTComputeFreeState();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}


int sicSTcomputePbInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTComputePb();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicSTupdateStateInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTupdateState();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicModelgetQInterface(char *fname)
{
  static int minrhs = 2, maxrhs = 2, minlhs = 1, maxlhs = 2;
  static int dimIndexDS1 = 1, dimIndexDS2 = 1, IndexDS;
  static int dimIndexVec1 = 1, dimIndexVec2 = 1, IndexVec;
  static int dimVal1 = 1, dimVal2 = 1, Val;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  /* Get IndexDS (1, int *)  */
  GetRhsVar(1, "i", &dimIndexDS1, &dimIndexDS2, &IndexDS);
  if (!(dimIndexDS1 * dimIndexDS2 > 0))
  {
    sciprint("Wrong parameter in ssicModelgetQ (number IndexDS has wrong size!)\r\n");
    Error(999);
    return 0;
  }
  /* Get IndexVec (2, int *)  */
  GetRhsVar(2, "i", &dimIndexVec1, &dimIndexVec2, &IndexVec);
  if (!(dimIndexVec1 * dimIndexVec2 > 0))
  {
    sciprint("Wrong parameter in ssicModelgetQ (number IndexVector has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "d", &dimVal1, &dimVal2, &Val);
  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicModelgetQ(stk(Val), *istk(IndexDS), *istk(IndexVec));
  /*  Return variable  */
  LhsVar(1) = 3;
  LhsVar(2) = 4;



  return 0;
}

int sicLagrangianLinearTIDSInterface(char *fname)
{
  static int minrhs = 8, maxrhs = 8, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nDof, Q0, Vel0, Mass, K, C, libname, fctname;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=8) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  /* Get nDof (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nDof);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (nDof has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get Q0 (2, double)  */
  GetRhsVar(2, "d", &dim1, &dim2, &Q0);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (Q0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get Vel0 (3, double)  */
  GetRhsVar(3, "d", &dim1, &dim2, &Vel0);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (Vel0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get Mass (4, double)  */
  GetRhsVar(4, "d", &dim1, &dim2, &Mass);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (Mass has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get K (5, double)  */
  GetRhsVar(5, "d", &dim1, &dim2, &K);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (K has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get C (6, double)  */
  GetRhsVar(6, "d", &dim1, &dim2, &C);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (C has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (7, char)  */
  GetRhsVar(7, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    printf("-->%d %d\n", dim1, dim2);
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (8, char)  */
  GetRhsVar(8, "c", &dim1, &dim2, &fctname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearTIDS (fctname has wrong size!)\r\n");
    Error(999);
    return 0;
  }


  CreateVar(9, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicLagrangianLinearTIDS(*istk(nDof), stk(Q0), stk(Vel0), stk(Mass), stk(K), stk(C), cstk(libname), cstk(fctname));
  /*  Return variable  */
  LhsVar(1) = 9;

  return 0;
}

int sicInteractionLLRInterface(char *fname)
{
  static int minrhs = 8, maxrhs = 8, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int name, nbDS, DS, nbRel, H, b, lawtype, e;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=8) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get name (1, char)  */
  GetRhsVar(1, "c", &dim1, &dim2, &name);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteractionLLR (name has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get nbDS (2, int)  */
  GetRhsVar(2, "i", &dim1, &dim2, &nbDS);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicInteractionLLR (nbDS has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get DS (3, int)  */
  GetRhsVar(3, "i", &dim1, &dim2, &DS);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteractionLLR (DS has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get nbRel (4, int)  */
  GetRhsVar(4, "i", &dim1, &dim2, &nbRel);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicInteractionLLR (nbRel has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get H (5, double)  */
  GetRhsVar(5, "d", &dim1, &dim2, &H);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteractionLLR (H has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get b (6, double)  */
  GetRhsVar(6, "d", &dim1, &dim2, &b);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteractionLLR (b has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get lawtype (7, char)  */
  GetRhsVar(7, "c", &dim1, &dim2, &lawtype);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteractionLLR (name has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get e (8, double)  */
  GetRhsVar(8, "d", &dim1, &dim2, &e);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicInteractionLLR (b has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(9, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicInteractionLLR(cstk(name), *istk(nbDS), istk(DS), *istk(nbRel), stk(H), stk(b), cstk(lawtype), *stk(e));

  /*  Return variable  */
  LhsVar(1) = 9;

  return 0;
}

int sicNSDSModelInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int isBVP, t0, T;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get isBVP (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &isBVP);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNSDSModel (isBVP has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get t0 (2, double)  */
  GetRhsVar(2, "d", &dim1, &dim2, &t0);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNSDSModel (t0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get T (3, double)  */
  GetRhsVar(3, "d", &dim1, &dim2, &T);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNSDSModel (T has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicNSDSModel(*istk(isBVP), *stk(t0), *stk(T));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicStrategyTimeSteppingInterface(char *fname)
{
  static int minrhs = 4, maxrhs = 4, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int h, theta, maxiter, tolerance;
  static int dimo1 = 1, dimo2 = 1, st;

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get  h (1, double)  */
  GetRhsVar(1, "d", &dim1, &dim2, &h);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicStrategyTimeStepping (h has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get  h (2, theta)  */
  GetRhsVar(2, "d", &dim1, &dim2, &theta);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicStrategyTimeStepping (theta has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get  h (3, maxiter)  */
  GetRhsVar(3, "d", &dim1, &dim2, &maxiter);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicStrategyTimeStepping (maxiter has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get  h (4, tolerance)  */
  GetRhsVar(4, "d", &dim1, &dim2, &tolerance);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicStrategyTimeStepping (tolerance has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(5, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicStrategyTimeStepping(*stk(h), stk(theta), *stk(maxiter), *stk(tolerance));

  /*  Return variable  */
  LhsVar(1) = 5;


  return 0;
}



/***************************************************
 * Gateways for Scilab
 ***************************************************/

int C2F(SiconosGateway)()
{
  gate_function function[] = {sicLoadModelInterface,
                              sicInitStrategyInterface,
                              sicTimeGetHInterface,
                              sicTimeGetNInterface,
                              sicTimeGetKInterface,
                              sicSTNextStepInterface,
                              sicSTComputeFreeStateInterface,
                              sicSTcomputePbInterface,
                              sicSTupdateStateInterface,
                              sicModelgetQInterface,
                              sicLagrangianLinearTIDSInterface,
                              sicInteractionLLRInterface,
                              sicNSDSModelInterface,
                              sicStrategyTimeSteppingInterface
                             };

  char *name[] = {"sicLoadModel",
                  "sicInitStrategy",
                  "sicTimeGetH",
                  "sicTimeGetN",
                  "sicTimeGetK",
                  "sicSTNextStep",
                  "sicSTComputeFreeState",
                  "sicSTcomputePb",
                  "sicSTupdateState",
                  "sicModelgetQ"
                  "sicLagrangianLinearTIDS",
                  "sicInteractionLLR",
                  "sicNSDSModel",
                  "sicStrategyTimeStepping"
                 };

  Rhs = Max(0, Rhs);
  sci_gateway(name[Fin - 1], function[Fin - 1]);
  return 0;
}

