#include "stack-c.h"
#include <string.h>
#include "machine.h"
#include "sun/link.h"
#include "robotmodel.h"

#define CreateVar1(n,ct,mx,nx,lx)\
 if(! C2F(createvar)((c_local=n,&c_local),ct,mx,nx,lx, 1L))\
        { return ;  }



/***************************************************
 * interfaces
 ***************************************************/

void modele_nddl(int *nddl)
{
  *nddl = NDDL;
}

int interface_nddl(char *fname)
{
  static int un = 1, nddl;

  /*    Define minls=1,  maxlhs   */
  static int minlhs = 1,  maxlhs = 1;

  /*   Check lhs   */
  CheckLhs(minlhs, maxlhs) ;

  CreateVar(1, "d", &un, &un, &nddl);

  modele_nddl(istk(nddl));

  LhsVar(1) = 1;
  return 0;
}

int interface_frottements(char *fname)
{
  static int un = 1, nddl = NDDL;
  static int n, nbis;
  static int q, qdot, F;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 2, maxlhs = 1, maxrhs = 2;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  GetRhsVar(2, "d", &n, &nbis, &qdot);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "d", &nddl, &un, &F);

  modele_frottements(stk(q), stk(qdot), stk(F));

  LhsVar(1) = 3;
  return 0;
}

int interface_coriolis(char *fname)
{
  static int nddl = NDDL;
  static int n, nbis;
  static int q, qdot, N;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 2, maxlhs = 1, maxrhs = 2;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  GetRhsVar(2, "d", &n, &nbis, &qdot);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(3, "d", &nddl, &nddl, &N);

  modele_coriolis(stk(q), stk(qdot), stk(N));

  LhsVar(1) = 3;
  return 0;
}

int interface_gravite(char *fname)
{
  static int un = 1, nddl = NDDL;
  static int n, nbis;
  static int q, G;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 1, maxlhs = 1, maxrhs = 1;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "d", &nddl, &un, &G);

  modele_gravite(stk(q), stk(G));

  LhsVar(1) = 2;
  return 0;
}


int interface_inertie(char *fname)
{
  static int nddl = NDDL;
  static int n, nbis;
  static int q, M;

  /*    Define minls=1, maxlhs, minrhs, maxrhs   */
  static int minlhs = 1, minrhs = 1, maxlhs = 1, maxrhs = 1;

  /*   Check rhs and lhs   */
  CheckRhs(minrhs, maxrhs) ;
  CheckLhs(minlhs, maxlhs) ;

  GetRhsVar(1, "d", &n, &nbis, &q);
  if (n * nbis != NDDL)
  {
    sciprint("Wrong size!\r\n");
    Error(999);
    return 0;
  }

  CreateVar(2, "d", &nddl, &nddl, &M);

  modele_inertie(stk(q), stk(M));

  LhsVar(1) = 2;
  return 0;
}



/***************************************************
 * Gateways for Scilab
 ***************************************************/


typedef int (*gate_function)(char *);

int C2F(gateway_lagrangien)()
{
  gate_function function[] = {interface_coriolis,
                              interface_gravite,
                              interface_inertie,
                              interface_frottements,
                              interface_nddl
                             };
  char *name[] = {"coriolis",
                  "gravite",
                  "inertie",
                  "frottements",
                  "nddl"
                 };

  Rhs = Max(0, Rhs);
  sci_gateway(name[Fin - 1], function[Fin - 1]);
  return 0;
}

