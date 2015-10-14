#include "stack-c.h"
#include <string.h>
#include "machine.h"
#include "sun/link.h"
#include "RobotModel.h"


/***************************************************
 * interfaces
 ***************************************************/

void Ndof(int *nddl)
{
  *nddl = NDOF;
}

int NdofInterface(char *fname)
{
  static int un = 1, nddl;

  /*    Define minls=1,  maxlhs   */
  static int minlhs = 1,  maxlhs = 1;

  /*   Check lhs   */
  CheckLhs(minlhs, maxlhs) ;

  CreateVar(1, "d", &un, &un, &nddl);

  Ndof(istk(nddl));

  LhsVar(1) = 1;
  return 0;
}

int FrictionInterface(char *fname)
{
  static int un = 1, nddl = NDOF;
  static int n, nbis;
  static int q, qdot, F;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 2, maxlhs = 1, maxrhs = 2;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDOF)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  GetRhsVar(2, "d", &n, &nbis, &qdot);
  if (n * nbis != NDOF)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "d", &nddl, &un, &F);

  Friction(stk(q), stk(qdot), stk(F));

  LhsVar(1) = 3;
  return 0;
}

int InertiaInterface(char *fname)
{
  static int ndof = NDOF;
  static int n, nbis;
  static int q, M;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 1, maxlhs = 1, maxrhs = 1;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDOF)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "d", &ndof, &ndof, &M);

  Inertia(stk(M), stk(q));

  LhsVar(1) = 2;
  return 0;
}

int NLEffectsInterface(char *fname)
{
  static int one = 1, ndof = NDOF;
  static int n, nbis;
  static int q, qdot, N;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 2, maxlhs = 1, maxrhs = 2;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDOF)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  GetRhsVar(2, "d", &n, &nbis, &qdot);
  if (n * nbis != NDOF)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "d", &ndof, &one, &N);

  NLEffects(stk(N), stk(q), stk(qdot));

  LhsVar(1) = 3;
  return 0;
}


/***************************************************
 * Gateways for Scilab
 ***************************************************/


typedef int (*gate_function)(char *);

int C2F(LagrangianGateway)()
{
  gate_function function[] = {NdofInterface,
                              FrictionInterface,
                              InertiaInterface,
                              NLEffectsInterface
                             };
  char *name[] = {"Ndof",
                  "Friction",
                  "Inertia",
                  "NLEffects"
                 };

  Rhs = Max(0, Rhs);
  sci_gateway(name[Fin - 1], function[Fin - 1]);
  return 0;
}

