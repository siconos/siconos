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
// @file ActuationModel/NoDynamics/TaskFunctionControl/ControlLaw.cpp
// @author Florence Billet
//
// Affiliation(s): INRIA, team BIPOP
//
// Email(s): Florence.Billet@inria.fr
//
// @brief Compute the voltage required to follow the reference trajectory
//
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <iostream>
#include <string>

#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/ublas/io.hpp>

extern "C" {
#include "Trajectory.h"
#include "TaskFunctionDefinition.h"

#include "LagrangianModel.h"
#include "SomeDefinitions.h"
}

#include "KernelSomeDefinitions.hpp"
#include "Utils.hpp"

//#include "../../../LagrangianDynamics/Complete/util.hpp" //a mettre ailleurs un jour

SICONOS_EXPORT void controlLaw(double * t, double * q, double * qdot, int * NDOF, int * NCONT, double * torques)
{
  int ndof;
  int contacts;
  char lagrangianModelName[20] = "";

  getLagrangianModelName(lagrangianModelName);
  getActiveDOFsize(&ndof);
  vector<int> activeDOF(ndof);
  if (ndof > 0)
    getActiveDOF(&(activeDOF[0]));
  else
    ndof = *NDOF;

  ///////////////////////////////////////////
  // Computation of the Lagrangian model
  ///////////////////////////////////////////


  matrix<double, column_major> M(*NDOF, *NDOF);
  vector<double> N(*NDOF);

  if (strcmp(lagrangianModelName, "Human36") == 0)
  {
    double L[31];
    double addl[84];
    double mass;
    GetModelMass(&mass);
    GetAnatomicalLengths(L);
    GetTag2JointLengths(addl);
    InertiaH36(&(M(0, 0)), q, L, addl, mass);
    NLEffectsH36(&(N[0]), q, qdot, L, addl, mass);
  }
  else
  {
    Inertia(&(M(0, 0)), q);
    NLEffects(&(N[0]), q, qdot);
  }

  //////////////////////////////////////////////////
  //Additionnal forces
  /////////////////////////////////////////////////

  vector<double> fAdditionnal(*NDOF);
  if (strcmp(lagrangianModelName, "PA10") == 0)
    Friction(&(fAdditionnal[0]), q, qdot);
  else if (strcmp(lagrangianModelName, "RX90") == 0)
    SpringForce(&(fAdditionnal[0]), q);
  else
    fAdditionnal.clear();

  ///////////////////////////////////////////////////
  // Limitation to the selected degrees of freedom
  ///////////////////////////////////////////////////
  reduceMatrix(M, activeDOF, activeDOF);
  N = reduceVector(N, activeDOF);
  fAdditionnal = reduceVector(fAdditionnal, activeDOF);

  //////////////////////////////////////////////////
  // Proportionnel-derivee in the task space
  //////////////////////////////////////////////////

  vector<double> sDesired(*NDOF);
  vector<double> sdotDesired(*NDOF);
  vector<double> sddotDesired(*NDOF);

  trajectory(t, &(sDesired[0]), &(sdotDesired[0]), &(sddotDesired[0]), &contacts);

  //v = Kp*(s_desiree-TaskFunction(q))+Kv*(sdot_desiree-H*qdot)-h+sddot_desiree;
  vector<double> TF(*NDOF);
  vector<double> h(*NDOF);
  matrix<double, column_major> H(*NDOF, *NDOF);

  TaskFunction(&(TF[0]), q);
  TaskJacobian(&(H(0, 0)), q);
  TaskNLEffects(&(h[0]), q, qdot);

  double Kp = 100.0;
  double Kv = 30.0;
  vector<double, array_adaptor<double> > tmp_cast(*NDOF, array_adaptor<double> (*NDOF, qdot));
  sddotDesired += Kp * (sDesired - TF) + Kv * (sdotDesired - prod(H, tmp_cast)) - h;

  /////////////////////////////////////////////////////////////
  // Generalized Torques to make the realisation of the command
  /////////////////////////////////////////////////////////////

  //qddot = inv(H)*v;
  vector<int> ipiv(*NDOF);
  boost::numeric::bindings::atlas::getrf(H, ipiv);
  boost::numeric::bindings::atlas::getri(H, ipiv);
  sddotDesired = prod(H, sddotDesired);

  //D = M*qddot(ActiveDOF) + N + FAdditionnal;
  sddotDesired = reduceVector(sddotDesired, activeDOF);

  N += prod(M, sddotDesired) + fAdditionnal;

  //nmot = sum(ActiveDOF<=size(q, 1)); en fait ici nmot = ndof
  //Torques = zeros(q);
  //Torques(ActiveDOF(1:nmot)) = D(1:nmot);
  for (int i = 0; i < *NDOF; i++)
    torques[i] = 0;
  expandVector(torques, N, activeDOF);

}
