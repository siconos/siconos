/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to the modeling, the simulation and the control
 * of non smooth dynamical systems
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include <iostream>
#include "Model.h"
#include "check.h"

#include "LagrangianDS.h"
#include "LinearDS.h"

#include "LagrangianLinearTIDS.h"
#include "LagrangianR.h"

#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"

#include "TimeStepping.h"
#include "Moreau.h"
#include "LCP.h"

#include <libxml/parser.h>
#include "SiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"
#include <math.h>
#include <stdio.h>
#include "LCP.h"

#include <sys/time.h>

#include "SiconosApiC.h"



using namespace std;



//
// TODO: add a structure with array and manage index
//


Model *GLOB_MODEL;
NonSmoothDynamicalSystem *GLOB_NSDS;
Strategy *GLOB_STRATEGIE;
TimeDiscretisation *GLOB_TIME;
vector<DynamicalSystem *> GLOB_VECTOR_DS;
vector<Interaction*> GLOB_VECTOR_INTERACTION;


extern "C" int sicLoadModel(char ModelXmlFile[])
{
  int ret = SIC_OK;

  try
  {

    GLOB_MODEL = new Model(ModelXmlFile);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLoadModel" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicInitStrategy()
{
  int ret = SIC_OK;

  try
  {
    if (GLOB_MODEL != NULL)
    {
      GLOB_STRATEGIE = GLOB_MODEL->getStrategyPtr();
      GLOB_STRATEGIE->initialize();
    }
    else
    {
      RuntimeException::selfThrow("siconos/C:: sicInitStrategy failed");
    }

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicInitStrategy" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicTimGetH(double *H)
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    *H = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getH();
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int  sicTimeGetN(int *N)
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    *N = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getNSteps();
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int  sicTimeGetK(int *K)
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    *K = GLOB_STRATEGIE->getTimeDiscretisationPtr()->getK();
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int sicSTNextStep()
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    GLOB_STRATEGIE->nextStep();
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int sicSTComputeFreeState()
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    GLOB_STRATEGIE->computeFreeState();
  else
    ret = SIC_ERROR;

  return ret;
}


extern "C" int sicSTComputePb()
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    GLOB_STRATEGIE->computeOneStepNSProblem();
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int sicSTnewtonSolve(double criterion, int maxIter)
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    GLOB_STRATEGIE->newtonSolve(criterion, maxIter);
  else
    ret = SIC_ERROR;

  return ret;
}

extern "C" int sicSTupdateState()
{
  int ret = SIC_OK;

  if (GLOB_STRATEGIE != NULL)
    GLOB_STRATEGIE->update();
  else
    ret = SIC_ERROR;

  return ret;
}


extern "C" void sicDebug(int *ret)
{
  // GLOB_STRATEGIE->update();
  cout << "--- Debug" << endl;

  // LagrangianDS* ball1 = static_cast<LagrangianDS*> (GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));

  //   vector<Interaction*> vIS= (GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getInteractions());

  //   for (int index = 0; index < vIS.size(); index++)
  //     {
  //       vIS[index]->display();
  //     }

  *ret = 0;
}

extern "C" int sicModelgetQ(double *value, int indexDS, int indexVector)
{

  int ret = SIC_OK;
  int size;

  try
  {

    LagrangianDS* system = static_cast<LagrangianDS*>(GLOB_MODEL->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(indexDS));

    if (system != NULL)
    {
      size = (system->getQ()).size();
      if ((indexVector >= size) || (indexVector < 0))
        RuntimeException::selfThrow("siconos/C:: sicModelgetQ failed");
      else
      {
        *value = system->getQ()(indexVector);
      }
    }

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicModelgetQ" << endl;
    ret = SIC_ERROR;
  }

  return ret;

}

extern "C" int sicLagrangianLinearTIDS(int nDof, double *Q0, double *Vel0, double *Mass, double *K, double *C, char *libname, char * fctname)
{
  int nId, i, j;

  try
  {
    // Id of LagrangianLinearTIDS
    nId = GLOB_VECTOR_DS.size();

    // TODO: parameters verification

    // Vectors and Matrix creation
    SimpleVector vQ0(nDof);
    SimpleVector vVel0(nDof);
    SiconosMatrix mMass(nDof, nDof);
    SiconosMatrix mK(nDof, nDof);
    SiconosMatrix mC(nDof, nDof);

    // Vectors and Matrix initialisation with function parameters
    // Is there a solution less stupid ?
    for (i = 0; i < nDof; i++)
    {
      vQ0.setValue(i, Q0[i]);
      vVel0.setValue(i, Vel0[i]);
      for (j = 0; j < nDof; j++)
      {
        mMass.setValue(i, j, Mass[j + i * nDof]);
        mK.setValue(i, j, Mass[j + i * nDof]);
        mC.setValue(i, j, Mass[j + i * nDof]);
      }
    }

    // Create the object
    DynamicalSystem *DS = new LagrangianLinearTIDS(nId, nDof, vQ0, vVel0, mMass, mK, mC);

    // Set the Fext
    static_cast<LagrangianDS*>(DS)->setComputeFExtFunction(libname, fctname);

    // Push the LagrangianLinearTIDS on global DS vectors
    GLOB_VECTOR_DS.push_back(DS);

    cout << "Create LagrangianLinearTIDS Id " << nId << endl;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLagrangianLinearTIDS" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

extern "C" int sicLagrangianDS(int nDof, double *Q0, double *Vel0)
{
  int nId, i;

  try
  {
    // Id of LagrangianDS
    nId = GLOB_VECTOR_DS.size();

    // Vectors and Matrix creation
    SimpleVector vQ0(nDof);
    SimpleVector vVel0(nDof);

    // Vectors initialisation with function parameters
    // Is there a solution less stupid ?
    for (i = 0; i < nDof; i++)
    {
      vQ0.setValue(i, Q0[i]);
      vVel0.setValue(i, Vel0[i]);
    }
    // Create the object
    LagrangianDS * DS = new LagrangianDS(nId, nDof, vQ0, vVel0);
    /*
    DS->setComputeFIntFunction(libname,fctFInt);
    DS->setComputeJacobianQFIntFunction(libname,fctJacFInt);
    DS->setComputeJacobianVelocityFIntFunction(libname, fctJacFInt);
    DS->setComputeFExtFunction(libname, fctFext);
    */
    // Push the LagrangianDS on global DS vectors
    GLOB_VECTOR_DS.push_back(DS);

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLagrangianLinearDS" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

extern "C" int sicSetComputeMassFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeMassFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeMassFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeNNLFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeNNLFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeNNLFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeJacobianQNNLFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQNNLFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeJacobianQNNLFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeJacobianQNNLFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int  sicSetComputeJacobianVelocityNNLFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)-> setComputeJacobianVelocityNNLFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in  sicSetComputeJacobianVelocityNNLFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeFIntFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeFIntFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeFIntFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeJacobianQFIntFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQFIntFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeJacobianQFIntFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeJacobianQFIntFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeJacobianVelocityFIntFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)->setComputeJacobianVelocityFIntFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeJacobianVelocityFIntFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int  sicSetComputeFExtFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    // DS verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    if (GLOB_VECTOR_DS[nIdDs] == NULL)
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = GLOB_VECTOR_DS[nIdDs];
    static_cast<LagrangianDS*>(DS)-> setComputeFExtFunction(libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in  setComputeFExtFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicInteraction(char *name, int nbDS, int *DS, int nbRel)
{
  int nId, i, dimDS = 0;

  try
  {

    // Id of Interaction
    nId = GLOB_VECTOR_INTERACTION.size();

    // TODO: parameters verification
    int  nDSMax = GLOB_VECTOR_DS.size();
    if (nDSMax < nbDS)
      RuntimeException::selfThrow("siconos/C:: sicInteraction failed");

    // Compute sum of DS size
    for (i = 0; i < nbDS; i++)
    {
      dimDS += static_cast<LagrangianDS *>(GLOB_VECTOR_DS[DS[i]])->getNdof();
    }

    vector <DynamicalSystem *> *dsConcerned = new  vector<DynamicalSystem *>(nbDS);
    for (i = 0; i < nbDS; i++)
    {
      (*dsConcerned)[i] = GLOB_VECTOR_DS[DS[i]];
    }

    // Create object
    Interaction *interaction = new  Interaction(name, nId, nbRel, dsConcerned);
    // interaction->display();

    // Push the interaction on global DS vectors
    GLOB_VECTOR_INTERACTION.push_back(interaction);

    cout << "Create Interaction Id " << nId << endl;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicInteraction" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}



extern "C" int sicLagrangianLinearR(int nIdInteraction, double *H, double *b)
{
  int nId, i, j, dimDS = 0, nbRel = 0;


  try
  {

    // Retrieve Interaction
    int nSize = GLOB_VECTOR_INTERACTION.size();
    if ((nIdInteraction < 0) || (nIdInteraction > nSize))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR failed");

    Interaction *interaction = GLOB_VECTOR_INTERACTION[nIdInteraction];

    // Compute sum of DS size
    dimDS = interaction->getSizeOfDS();
    nbRel = interaction->getNInteraction();

    if ((dimDS <= 0) || (nbRel <= 0))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR  failed");

    // Vectors and Matrix creation
    SiconosMatrix mH(nbRel, dimDS);
    SimpleVector  vb(nbRel);

    // Vectors and Matrix initialisation with function parameters
    // Is there a solution less stupid ?
    vb.setValue(0, 0);
    for (i = 0; i < nbRel; i++)
    {
      vb(i) = b[i];
      for (j = 0; j < dimDS; j++)
      {
        mH.setValue(i, j, H[i * nbRel + j]);
      }
    }
    cout << "H(" << dimDS << ")=" << mH << endl;

    Relation *relation = new LagrangianLinearR(mH, vb);
    interaction->setRelationPtr(relation);

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLagrangianLinearR" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

extern "C" int sicLagrangianR(int nIdInteraction, char *relationType, char *funcH, char *funcG)
{
  int nId = SIC_OK;


  try
  {
    // Retrieve Interaction
    int nSize = GLOB_VECTOR_INTERACTION.size();
    if ((nIdInteraction < 0) || (nIdInteraction > nSize))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianR failed");

    Interaction *interaction = GLOB_VECTOR_INTERACTION[nIdInteraction];

    if (interaction == NULL)
      RuntimeException::selfThrow("siconos/C:: sicLagrangianR failed");

    // Why a vector ? See Franck or Vincent...
    vector<string> vG;
    vG.reserve(1);
    vG.push_back(funcG);

    Relation * relation = new LagrangianR(relationType, funcH, vG, interaction);
    interaction->setRelationPtr(relation);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLagrangianLinearR" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}


extern "C" int sicNewtonImpactLawNSL(int nIdInteraction, double e)
{
  int nId = SIC_OK;

  try
  {

    // Retrieve Interaction
    int nSize = GLOB_VECTOR_INTERACTION.size();
    if ((nIdInteraction < 0) || (nIdInteraction > nSize))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR failed");

    Interaction *interaction = GLOB_VECTOR_INTERACTION[nIdInteraction];

    NonSmoothLaw *nslaw = new NewtonImpactLawNSL(e);

    interaction->setNonSmoothLawPtr(nslaw);

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in ssicNewtonImpactLawNSL" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}


extern "C" int sicNonSmoothDynamicalSystem(int isBVP)
{
  int nId = SIC_OK;

  try
  {

    if ((isBVP < 0) || (isBVP > 1))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad isBVP");
    // Create NSDS system connected with DS and Interactions
    GLOB_NSDS = new NonSmoothDynamicalSystem(isBVP);

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicNonSmoothDynamicalSystem" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

extern "C" int sicModel(double t0, double T)
{
  int nId = SIC_OK;

  try
  {

    // Parameters verification
    if ((t0 < 0.0) || (T < 0.0))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad time parameters");

    if (GLOB_VECTOR_DS.size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to DS empty");


    if (GLOB_VECTOR_INTERACTION.size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to Interaction empty");
    if (GLOB_NSDS == NULL)
      RuntimeException::selfThrow("siconos/C:: sicModel failed due to NSDS empty");

    GLOB_NSDS->setDynamicalSystems(GLOB_VECTOR_DS);
    GLOB_NSDS->setInteractions(GLOB_VECTOR_INTERACTION);

    // Create the model connected with NSDS
    GLOB_MODEL = new Model(t0, T);
    GLOB_MODEL->setNonSmoothDynamicalSystemPtr(GLOB_NSDS);
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicInteraction" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

extern "C" int sicStrategyTimeStepping(double h)
{
  int nId = SIC_OK;

  try
  {

    GLOB_STRATEGIE = new TimeStepping(GLOB_MODEL);

    // Time discretisation
    GLOB_TIME = new  TimeDiscretisation(h, GLOB_STRATEGIE);

    GLOB_STRATEGIE->initialize();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicStrategyTimeStepping" << endl;
    nId = SIC_ERROR;
  }

  return nId;

}


extern "C" int sicOneStepIntegratorMoreau(double *theta)
{
  int nId = SIC_OK, i;

  try
  {
    int  dsNumber = GLOB_VECTOR_DS.size();

    if (dsNumber == 0)
      RuntimeException::selfThrow("siconos/C:: ssicStrategyTimeStepping failed due to DS empty");

    // One Step integrator s
    vector <OneStepIntegrator *> vOSI;
    vOSI.resize(dsNumber, NULL);
    TimeDiscretisation *t = GLOB_STRATEGIE->getTimeDiscretisationPtr();
    for (i = 0; i < dsNumber; i++)
      vOSI[i] = new Moreau(t, GLOB_VECTOR_DS[i], theta[i]);

    GLOB_STRATEGIE->setOneStepIntegrators(vOSI);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicOneStepIntegratorMoreau" << endl;
    nId = SIC_ERROR;
  }

  return nId;

}

extern "C" int sicOneStepNSProblemLCP(char *solverName, int maxiter, double tolerance)
{
  int nId = SIC_OK;

  try
  {
    // One Step problem
    OneStepNSProblem * onsspb = new LCP(GLOB_STRATEGIE, solverName, maxiter, tolerance);
    GLOB_STRATEGIE->setOneStepNSProblemPtr(onsspb);

    cout << "=== End of model loading === " << endl;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicOneStepNSProblemLCP" << endl;
    nId = SIC_ERROR;
  }

  return nId;

}

extern "C" int sicClean()
{
  unsigned int i;

  delete GLOB_MODEL;
  delete GLOB_NSDS;
  delete GLOB_STRATEGIE;

  for (i = 0; i < GLOB_VECTOR_DS.size(); i++)
  {
    delete GLOB_VECTOR_DS[i];
  }

  for (i = 0; i < GLOB_VECTOR_INTERACTION.size(); i++)
  {
    delete GLOB_VECTOR_INTERACTION[i];
  }

  return SIC_OK;
}
