#include "SiconosSci.h"


int sicLoadModelInterface(char *fname)
{
  static int minrhs = 1, maxrhs = 1, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int dimo1 = 1, dimo2 = 1, st;
  static int ModelXmlFile;

#ifdef _DEBUG
  printf("sicLoadModelInterface\n");
#endif

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

int sicInitSimulationInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicInitSimulationInterface\n");
#endif


  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicInitSimulation();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicTimeGetHInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 2;
  static int dim1 = 1, dim2 = 1, H;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicTimeGetHInterface\n");
#endif


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

#ifdef _DEBUG
  printf("sicTimeGetNInterface\n");
#endif


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

int sicSTNextStepInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTNextStepInterface\n");
#endif

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

int sicSTAdvanceToEventInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTNextStepInterface\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTAdvanceToEvent();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}



int sicSTSaveInMemoryInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTSaveInMemory\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTSaveInMemory();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicSTComputeOneStepInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTAdvanceToEventInterface\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicSTComputeOneStep();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

int sicSTnewtonSolveInterface(char *fname)
{
  int sicSTnewtonSolve(double criterion, int maxIter);
  static int minrhs = 2, maxrhs = 2, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int criterion, maxIter;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTnewtonSolveInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get criterion (1, double)  */
  GetRhsVar(1, "d", &dim1, &dim2, &criterion);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in  sicSTnewtonSolve (criterion has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get maxIter (2, int)  */
  GetRhsVar(2, "i", &dim1, &dim2, &maxIter);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in  sicSTnewtonSolve (maxIter has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSTnewtonSolve(*stk(criterion), *istk(maxIter));

  /*  Return variable  */
  LhsVar(1) = 3;

  return 0;
}


int sicSTupdateStateInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSTupdateStateInterface\n");
#endif

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

#ifdef _DEBUG
  printf("sicModelgetQInterface\n");
#endif


  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  /* Get IndexDS (1, int *)  */
  GetRhsVar(1, "i", &dimIndexDS1, &dimIndexDS2, &IndexDS);
  if (!(dimIndexDS1 * dimIndexDS2 > 0))
  {
    sciprint("Wrong parameter in sicModelgetQ (number IndexDS has wrong size!)\r\n");
    Error(999);
    return 0;
  }
  /* Get IndexVec (2, int *)  */
  GetRhsVar(2, "i", &dimIndexVec1, &dimIndexVec2, &IndexVec);
  if (!(dimIndexVec1 * dimIndexVec2 > 0))
  {
    sciprint("Wrong parameter in sicModelgetQ (number IndexVector has wrong size!)\r\n");
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


#ifdef _DEBUG
  printf("sicLagrangianLinearTIDSInterface\n");
#endif

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

int sicLagrangianDSInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nDof, Q0, Vel0;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicLagrangianDSInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nDof);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicLagrangianDS (nDof has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get Q0 (2, double vector)  */
  GetRhsVar(2, "d", &dim1, &dim2, &Q0);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianDS (Q0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get Vel0 (3, double vector)  */
  GetRhsVar(3, "d", &dim1, &dim2, &Vel0);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianDS (Vel0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicLagrangianDS(*istk(nDof), stk(Q0), stk(Vel0));

  /*  Return variable  */
  LhsVar(1) = 4;


  return 0;
}



int sicSetMassInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetMassInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetMass (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetMass ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetMass (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeMassFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicSetNNLInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetNNLInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetNNL (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in ssicSetNNL ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetNNL (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeNNLFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicSetJacQNNLInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetJacQNNLInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetJacQNNL (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in ssicSetJacQNNL ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacQNNL (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeJacobianQNNLFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicSetJacVelNNLInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf(" sicSetJacVelNNLInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetJacVelNNL (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacVelNNL ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacVelNNL (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) =  sicSetComputeJacobianVelocityNNLFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicSetFIntInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetFIntInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetFInt (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetFInt (libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetFInt (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeFIntFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}

int sicSetJacQFIntInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetJacQFIntInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetJacQFInt (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacQFInt ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacQFInt (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) =  sicSetComputeJacobianQFIntFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}


int sicSetJacVelFIntInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetJacVelFIntInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetJacVelFInt (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacVelFInt ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetJacVelFInt (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeJacobianVelocityFIntFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}
int sicSetFExtInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, libname, func;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetFExtInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSetFExt  (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &libname);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in ssicSetFExt ( libname has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &func);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicSetFExt (func has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSetComputeFExtFunction(*istk(nId), cstk(libname), cstk(func));

  /*  Return variable  */
  LhsVar(1) = 4;

  return 0;
}


int sicInteractionInterface(char *fname)
{
  static int minrhs = 4, maxrhs = 4, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int name, nbDS, DS, nbRel;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicInteractionInterface\n");
#endif

  /* Check number of inputs (rhs=4) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get name (1, char)  */
  GetRhsVar(1, "c", &dim1, &dim2, &name);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteraction (name has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get nbDS (2, int)  */
  GetRhsVar(2, "i", &dim1, &dim2, &nbDS);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicInteraction (nbDS has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get DS (3, int)  */
  GetRhsVar(3, "i", &dim1, &dim2, &DS);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicInteraction (DS has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get nbRel (4, int)  */
  GetRhsVar(4, "i", &dim1, &dim2, &nbRel);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicInteraction (nbRel has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(5, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicInteraction(cstk(name), *istk(nbDS), istk(DS), *istk(nbRel));

  /*  Return variable  */
  LhsVar(1) = 5;

  return 0;
}


int sicLagrangianLinearRInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nIdInteraction, H, b;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicLagrangianLinearRInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nIdInteraction);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicLagrangianLinearR (nIdInteraction has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get H (2, double vector)  */
  GetRhsVar(2, "d", &dim1, &dim2, &H);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearR (H has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get b (3, double vector)  */
  GetRhsVar(3, "d", &dim1, &dim2, &b);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianLinearR (b has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  // OBSOLETE 2/3/07
  // *istk(st)= sicLagrangianLinearR(*istk(nIdInteraction),stk(H),stk(b));

  /*  Return variable  */
  LhsVar(1) = 4;


  return 0;
}

int sicLagrangianRInterface(char *fname)
{
  static int minrhs = 4, maxrhs = 4, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nId, relationType, funcH, funcG;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSetJacQNNLInterface\n");
#endif

  /* Check number of inputs (rhs=4) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nId);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicLagrangianR (nId has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get libname (2, char* )  */
  GetRhsVar(2, "c", &dim1, &dim2, &relationType);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianR (relationType has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get func (3, char* )  */
  GetRhsVar(3, "c", &dim1, &dim2, &funcH);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianR (funcH has wrong size!)\r\n");
    Error(999);
    return 0;
  }


  /* Get func (4, char* )  */
  GetRhsVar(4, "c", &dim1, &dim2, &funcG);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicLagrangianR (funcG has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(5, "i", &dimo1, &dimo2, &st);

  /* Call function */
  // OBSOLETE 2/3/07
  //  *istk(st)= sicLagrangianR(*istk(nId),cstk(relationType),cstk(funcH),cstk(funcG));

  /*  Return variable  */
  LhsVar(1) = 5;

  return 0;
}

int sicNewtonImpactLawNSLInterface(char *fname)
{
  static int minrhs = 2, maxrhs = 2, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int nIdInteraction, e;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicNewtonImpactNSL\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get nIdInteraction (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &nIdInteraction);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNewtonImpacNSL (nIdInteraction has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get e (2, double)  */
  GetRhsVar(2, "d", &dim1, &dim2, &e);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNewtonImpactNSL has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicNewtonImpactNSL(*istk(nIdInteraction), *stk(e));

  /*  Return variable  */
  LhsVar(1) = 3;

  return 0;
}

int sicNonSmoothDynamicalSystemInterface(char *fname)
{
  static int minrhs = 1, maxrhs = 1, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int isBVP;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicNonSmoothDynamicalSystemInterface\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get isBVP (1, int)  */
  GetRhsVar(1, "i", &dim1, &dim2, &isBVP);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicNonSmoothDynamicalSystem (isBVP has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicNonSmoothDynamicalSystem(*istk(isBVP));

  /*  Return variable  */
  LhsVar(1) = 2;

  return 0;
}


int sicModelInterface(char *fname)
{
  static int minrhs = 2, maxrhs = 2, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int isBVP, t0, T;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicModelInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  /* Get t0 (1, double)  */
  GetRhsVar(1, "d", &dim1, &dim2, &t0);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicModel (t0 has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get T (2, double)  */
  GetRhsVar(2, "d", &dim1, &dim2, &T);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicModel (T has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicModel(*stk(t0), *stk(T));

  /*  Return variable  */
  LhsVar(1) = 3;

  return 0;
}

int sicSimulationTimeSteppingInterface(char *fname)
{
  static int minrhs = 1, maxrhs = 1, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int h, theta, maxiter, tolerance;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicSimulationTimeSteppingInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get  h (1, double)  */
  GetRhsVar(1, "d", &dim1, &dim2, &h);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicSimulationTimeStepping (h has wrong size!)\r\n");
    Error(999);
    return 0;
  }


  CreateVar(2, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicSimulationTimeStepping(*stk(h));

  /*  Return variable  */
  LhsVar(1) = 2;

  return 0;
}

int sicOneStepIntegratorMoreauInterface(char *fname)
{
  static int minrhs = 1, maxrhs = 1, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int theta;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicOneStepIntegratorMoreauInterface\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get  theta (2, double vector)  */
  GetRhsVar(1, "d", &dim1, &dim2, &theta);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in ssicOneStepIntegratorMoreau (theta has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicOneStepIntegratorMoreau(stk(theta));

  /*  Return variable  */
  LhsVar(1) = 2;

  return 0;
}

int sicOneStepNSProblemLCPInterface(char *fname)
{
  static int minrhs = 3, maxrhs = 3, minlhs = 1, maxlhs = 1;
  static int dim1, dim2;
  static int solverName, maxiter, tolerance;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicOneStepNSProblemLCPInterface\n");
#endif

  /* Check number of inputs (rhs=3) and outputs (lhs=1) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  /* Get  solverName (1, char *)  */
  GetRhsVar(1, "c", &dim1, &dim2, &solverName);
  if (!(dim1 * dim2 > 0))
  {
    sciprint("Wrong parameter in sicOneStepNSProblemLCP (solverName has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get  maxiter (2, int)  */
  GetRhsVar(2, "i", &dim1, &dim2, &maxiter);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicOneStepNSProblemLCP (maxiter has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  /* Get  tolerance (3, double)  */
  GetRhsVar(3, "d", &dim1, &dim2, &tolerance);
  if (!(dim1 * dim2 == 1))
  {
    sciprint("Wrong parameter in sicOneStepNSProblemLCP (tolerance has wrong size!)\r\n");
    Error(999);
    return 0;
  }

  CreateVar(4, "i", &dimo1, &dimo2, &st);

  /* Call function */
  *istk(st) = sicOneStepNSProblemLCP(cstk(solverName), *istk(maxiter), *stk(tolerance));

  /*  Return variable  */
  LhsVar(1) = 4;


  return 0;
}

int sicCleanInterface(char *fname)
{
  static int minrhs = 0, maxrhs = 0, minlhs = 1, maxlhs = 1;
  static int dimo1 = 1, dimo2 = 1, st;

#ifdef _DEBUG
  printf("sicCleanInterface\n");
#endif

  /* Check number of inputs (rhs=1) and outputs (lhs=0) */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;


  CreateVar(1, "i", &dimo1, &dimo2, &st);
  /* Call function */
  *istk(st) = sicClean();
  /*  Return variable  */
  LhsVar(1) = 1;

  return 0;
}

