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
#include "../LagrangianModel.h"


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

extern int TagsInterface(char *fname);


/***************************************************
 * interfaces for Lagrangian Model
 ***************************************************/


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
  gate_function function[] = {TagsInterface
                             };
  char *name[] = {"Tags"
                 };

  Rhs = Max(0, Rhs);
  sci_gateway(name[Fin - 1], function[Fin - 1]);
  return 0;
}

