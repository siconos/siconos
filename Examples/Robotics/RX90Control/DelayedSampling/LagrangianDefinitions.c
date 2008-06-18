// Copyright (C) INRIA 1999-2008
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 2 as published
// by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//%
// @file LagrangianModel/RX90/SomeDefinitions.c
// @author Pierre-Brice Wieber
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Pierre-Brice.Wieber@inria.fr
//
// @brief Define and initialize the global variables MAXQ, MINQ,
// LagrangianModelName, ContactNames and ContactSolids
//

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "LagrangianDefinitions.h"

#define NDOF 6
#define NCONT 2
#define deg2rad M_PI/180.0


//Defintion of Bip joints limits

/**
 * Global variable corresponding to the max joints limit
 *
 * \var MAXQ (double) <em> dim NDOF </em>
 */
double MAXQ[6] = {167 * deg2rad, 137.5 * deg2rad, 142.5 * deg2rad, 270 * deg2rad, 112.5 * deg2rad, 270 * deg2rad};

void setMaxq(double * maxq)
{
  int i;

  for (i = 0; i < NDOF; i++)
    MAXQ[i] = maxq[i];
}

void getMaxq(double *maxq)
{
  int i;

  for (i = 0; i < NDOF; i++)
    maxq[i] = MAXQ[i];
}

/**
 * Global variable corresponding to the min joints limit
 *
 * \var MINQ (double) <em> dim NDOF </em>
 */
double MINQ[6] = { -167 * deg2rad, -137.5 * deg2rad, -142.5 * deg2rad, -270 * deg2rad, -112.5 * deg2rad, -270 * deg2rad};

void setMinq(double * minq)
{
  int i;

  for (i = 0; i < NDOF; i++)
    MINQ[i] = minq[i];
}

void getMinq(double *minq)
{
  int i;

  for (i = 0; i < NDOF; i++)
    minq[i] = MINQ[i];
}

//Definitions for Lagrangian Model

/**
 * Global variable storing the name of the lagrangian model
 *
 * \var LARANGIANMODELNAME (string)
 */
char *LAGRANGIANMODELNAME = "RX90";

void setLagrangianModelName(char* lagrangianModelName)
{
  LAGRANGIANMODELNAME = strdup(lagrangianModelName);
}

void getLagrangianModelName(char *lagrangianModelName)
{
  strcpy(lagrangianModelName, LAGRANGIANMODELNAME);
}

/**
 * Global variable storing the name of the contact points
 *
 * \var CONTACTNAMES (string) <em> dim NCONT </em>
 */
char *CONTACTNAMES[2] = {"", ""};

void setContactNames(char* contactNames, int *index)
{
  CONTACTNAMES[*index]  = strdup(contactNames);
}


void getContactNames(char *contactNames, int *index)
{
  strcpy(contactNames, CONTACTNAMES[*index]);
}

/**
 * Global variable storing the contact solids
 *
 * \var CONTACTSOLIDS (int) <em> dim 2*nbContactSolid </em>
 */
int CONTACTSOLIDS[2] = {1, 2};

void setContactSolids(int* contactSolids)
{
  int i;

  for (i = 0; i < 2; i++)
    CONTACTSOLIDS[i] = contactSolids[i];
}

void getContactSolids(int *contactSolids)
{
  int i;

  for (i = 0; i < 1; i++)
    contactSolids[i] = CONTACTSOLIDS[i];
}

void getNbContactSolids(int *nbSolids)
{
  *nbSolids = 1;
}


void SetModelSize(double *subjectSize) {}

void GetModelSize(double *subjectSize) {}

void SetModelMass(double *subjectMass) {}

void GetModelMass(double *subjectMass) {}

void SetAnatomicalLengths(double *anatLengths) {}

void GetAnatomicalLengths(double *anatLengths) {}

void SetTag2JointLengths(double *tagLengths) {}

void GetTag2JointLengths(double *tagLengths) {}

void Friction(double *F, double *q, double *qdot) {}

void TagsH36(double *T, double *q, double *L, double *addl, double mass) {}

void InertiaH36(double *M, double *q, double *L, double *addl, double mass) {}

void NLEffectsH36(double *N, double *q, double *qdot, double *L, double *addl, double mass) {}

void ContactHessianH36(double *H, double *q, double *qdot, double *L, double *addl) {}

void ContactJacobianH36(double *CJ, double *q, double *L, double *addl) {}

void ContactH36(double *CC, double *q, double *L, double *addl) {}

