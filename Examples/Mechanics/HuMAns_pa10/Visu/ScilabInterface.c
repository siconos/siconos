/* Copyright (C) INRIA 1999-2005
**
** This program is free software; you can redistribute it and/or modify it
** under the terms of the GNU General Public License version 2 as published
** by the Free Software Foundation.
**
** This program is distributed in the hope that it will be useful, but
** WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
** Public License for more details.
**
** You should have received a copy of the GNU General Public License along
** with this program; if not, write to the Free Software Foundation, Inc.,
** 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**
*/
/**
** @file LagrangianModel.h
** @author: Pierre-Brice Wieber
** Affiliation(s): INRIA, team BIPOP
** Email(s): Pierre-Brice.Wieber@inria.fr
**
** @brief Compute the complete Dynamics
**
**
*/


#include "stack-c.h"
#include <string.h>
#include "machine.h"
#include "os_specific/link.h"
#include "LagrangianModel.h"

#include <stdio.h>

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif


/***************************************************
 * Declarations for Lagrangian gateway             *
 ***************************************************/

extern int C2F(LagrangianGateway)();

/***************************************************
 * Declarations for Lagrangian Model               *
 ***************************************************/

extern int ContactInterface(char *fname);

extern int ContactHessianInterface(char *fname);

extern int ContactJacobianInterface(char *fname);

extern int InertiaInterface(char *fname);

extern int NLEffectsInterface(char *fname);

extern int JacobianNLEffectsInterface(char *fname);

extern int JacobianVelocityNLEffectsInterface(char *fname);

extern int TagsInterface(char *fname);



/***************************************************
 * interfaces for Lagrangian Model
 ***************************************************/

int ContactInterface(char *fname)
{
  static int one = 1, ncont = 3 * NCONT;
  static int n, nbis;
  static int q, CC;

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

  CreateVar(2, "d", &ncont, &one, &CC);

  Contact(stk(CC), stk(q));

  LhsVar(1) = 2;
  return 0;
}

int ContactHessianInterface(char *fname)
{
  static int one = 1, ncont = 3 * NCONT;
  static int n, nbis;
  static int q, qdot, H;

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

  CreateVar(3, "d", &ncont, &one, &H);

  ContactHessian(stk(H), stk(q), stk(qdot));

  LhsVar(1) = 3;
  return 0;
}

int ContactJacobianInterface(char *fname)
{
  static int ndof = NDOF, ncont = 3 * NCONT;
  static int n, nbis;
  static int q, J;

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

  CreateVar(2, "d", &ncont, &ndof, &J);

  ContactJacobian(stk(J), stk(q));

  LhsVar(1) = 2;
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



int JacobianNLEffectsInterface(char *fname)
{
  static int ndof = NDOF;
  static int n, nbis;
  static int q, qdot, NJ;

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

  CreateVar(3, "d", &ndof, &ndof, &NJ);

  JacobianNLEffects(stk(NJ), stk(q), stk(qdot));

  LhsVar(1) = 3;
  return 0;
}


int JacobianVelocityNLEffectsInterface(char *fname)
{
  static int ndof = NDOF;
  static int n, nbis;
  static int q, qdot, NVJ;

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

  CreateVar(3, "d", &ndof, &ndof, &NVJ);

  JacobianVelocityNLEffects(stk(NVJ), stk(q), stk(qdot));

  LhsVar(1) = 3;
  return 0;
}



int TagsInterface(char *fname)
{
  static int three = 3, ntags = NTAGS;
  static int n, nbis;
  static int q, T;

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

  CreateVar(2, "d", &ntags, &three, &T);

  Tags(stk(T), stk(q));

  LhsVar(1) = 2;
  return 0;
}




/***************************************************
 * Gateways for Scilab
 ***************************************************/

typedef int (*gate_function)(char *);

extern int sci_gateway(char *name, gate_function f);

int C2F(LagrangianGateway)()
{
  gate_function function[] = {ContactInterface,
                              ContactHessianInterface,
                              ContactJacobianInterface,
                              InertiaInterface,
                              NLEffectsInterface,
                              JacobianNLEffectsInterface,
                              JacobianVelocityNLEffectsInterface,
                              TagsInterface
                             };
  char *name[] = {"Contact",
                  "ContactHessian",
                  "ContactJacobian",
                  "Inertia",
                  "NLEffects",
                  "JacobianNLEffects",
                  "JacobianVelocityNLEffects",
                  "Tags"
                 };

  Rhs = Max(0, Rhs);
  sci_gateway(name[Fin - 1], function[Fin - 1]);
  return 0;
}

