/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#include "SiconosKernel.h"
#include "SiconosDataC.h"
#include "SiconosApiC.h"



//
// TODO: add a structure with array and manage index
//


//Model *GLOB_MODEL;
//NonSmoothDynamicalSystem *GLOB_NSDS;
//Simulation *GLOB_SIMULATION;
//TimeDiscretisation *GLOB_TIME;
//DynamicalSystemsSet *GLOB_SET_DS;
//InteractionsSet GLOB_SET_INTERACTION;
//CheckInsertDS GLOB_CHECK_DS;

DataC GLOB_DATA;

extern "C" int sicLoadModel(char ModelXmlFile[])
{
  int ret = SIC_OK;

  try
  {
    Model *ptrModel = new Model(ModelXmlFile) ;
    GLOB_DATA.setModelPtr(ptrModel);
    GLOB_DATA.setSimulationPtr(ptrModel->getSimulationPtr());
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

extern "C" int sicInitSimulation()
{
  int ret = SIC_OK;

  try
  {

    Model * prtModel = GLOB_DATA.getModelPtr();
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();



    prtSimul = prtModel->getSimulationPtr();
    prtSimul->initialize();

    SetDSPtr = new DynamicalSystemsSet();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicInitSimulation" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicTimeGetH(double *H)
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    *H = prtSimul->getTimeDiscretisationPtr()->getH();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicTimeGetH" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int  sicTimeGetN(int *N)
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    *N = prtSimul->getTimeDiscretisationPtr()->getNSteps();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicTimeGetNl" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicSTNextStep()
{
  int ret = SIC_OK;
  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    //RPG OUPS !!!
    ((TimeStepping*) prtSimul)->nextStep();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTNextStep" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicSTAdvanceToEvent()
{
  int ret = SIC_OK;
  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    //RPG OUPS !!!
    ((TimeStepping*) prtSimul)->advanceToEvent();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught insicSTAdvanceToEvent" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicSTSaveInMemory()
{
  int ret = SIC_OK;
  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    //RPG OUPS !!!
    prtSimul->saveInMemory();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTSaveInMemory" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicSTComputeOneStep()
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    ((TimeStepping*) prtSimul)->computeOneStep();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTComputeOneStep" << endl;
    ret = SIC_ERROR;
  }


  return ret;
}

extern "C" int sicSTnewtonSolve(double criterion, int maxIter)
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    //RPG OUPS !!!
    ((TimeStepping*) prtSimul)->newtonSolve(criterion, maxIter);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTnewtonSolve" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}


extern "C" int sicSTupdateState()
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    prtSimul->update(1); // WARNING FP: update needs an input, ie the derivative order used for lambda to compute R
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSTupdateState" << endl;
    ret = SIC_ERROR;
  }


  return ret;
}


extern "C" void sicDebug(int *ret)
{
  // GLOB_SIMULATION->update();
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

    Model * prtModel = GLOB_DATA.getModelPtr();

    LagrangianDS* system = static_cast<LagrangianDS*>(prtModel->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(indexDS));

    if (system != NULL)
    {
      size = (system->getQ()).size();
      if ((indexVector >= size) || (indexVector < 0))
        RuntimeException::selfThrow("siconos/C:: sicModelgetQ failed");
      else
      {
        printf("DEBUG --> %d %lf \n", indexVector, system->getQ()(indexVector));
        *value = system->getQ()(indexVector);
      }
    }
    else
    {
      RuntimeException::selfThrow("siconos/C:: sicModelgetQ failed");
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

    DynamicalSystemsSet* SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();
    CheckInsertDS* CheckDSPtr = GLOB_DATA.getCheckInsertDStPtr();

    // Id of LagrangianLinearTIDS
    nId = SetDSPtr->size();

    // TODO: parameters verification

    // Vectors and Matrix creation
    SimpleVector vQ0(nDof);
    SimpleVector vVel0(nDof);
    SimpleMatrix  mMass(nDof, nDof);
    SimpleMatrix  mK(nDof, nDof);
    SimpleMatrix mC(nDof, nDof);

    // Vectors and Matrix initialisation with function parameters
    // Is there a solution less stupid ?
    for (i = 0; i < nDof; i++)
    {
      vQ0(i) = Q0[i];
      vVel0(i) = Vel0[i];
      for (j = 0; j < nDof; j++)
      {
        mMass(i, j) = Mass[j + i * nDof];
        mK(i, j) = Mass[j + i * nDof];
        mC(i, j) = Mass[j + i * nDof];
      }
    }

    // Create the object
    //DynamicalSystem *DS = new LagrangianLinearTIDS(nId,nDof,vQ0,vVel0,mMass,mK,mC);

    // Push the LagrangianLinearTIDS on global DS vectors
    *CheckDSPtr = SetDSPtr->insert(new LagrangianLinearTIDS(nId, vQ0, vVel0, mMass, mK, mC));

    // Set the Fext
    (static_cast<LagrangianDS*>(*(CheckDSPtr->first)))->setComputeFExtFunction(libname, fctname);

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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // Id of LagrangianDS
    nId = SetDSPtr->size();

    // Vectors and Matrix creation
    SimpleVector vQ0(nDof);
    SimpleVector vVel0(nDof);

    // Vectors initialisation with function parameters
    // Is there a solution less stupid ?
    for (i = 0; i < nDof; i++)
    {
      vQ0(i) = Q0[i];
      vVel0(i) = Vel0[i];
    }
    // Create the object
    LagrangianDS * DS = new LagrangianDS(nId, vQ0, vVel0);
    /*
      DS->setComputeFIntFunction(libname,fctFInt);
      DS->setComputeJacobianQFIntFunction(libname,fctJacFInt);
      DS->setComputeJacobianVelocityFIntFunction(libname, fctJacFInt);
      DS->setComputeFExtFunction(libname, fctFext);
    */
    // Push the LagrangianDS on global DS vectors
    SetDSPtr->insert(DS);

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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }

    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQNNLFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)->setComputeJacobianNNLFunction(0, libname, func);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)-> setComputeJacobianNNLFunction(1, libname, func);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQFIntFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianQFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)->setComputeJacobianFIntFunction(0, libname, func);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)->setComputeJacobianFIntFunction(1, libname, func);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    if (SetDSPtr->getDynamicalSystemPtr(nIdDs) == NULL)
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getDynamicalSystemPtr(nIdDs);
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
  int nId, dimDS = 0;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();
    InteractionsSet * InteractionSetPtr = GLOB_DATA.getInteractionsSetPtr();

    // Id of Interaction
    nId = InteractionSetPtr->size();

    // TODO: parameters verification
    int  nDSMax = SetDSPtr->size();
    if (nDSMax < nbDS)
      RuntimeException::selfThrow("siconos/C:: sicInteraction failed");

    // Compute sum of DS size and add DS into the set
    DSIterator it;
    DynamicalSystemsSet dsConcerned;
    for (it = SetDSPtr->begin(); it != SetDSPtr->end(); ++it)
    {
      dimDS += (*it)->getDim();
      dsConcerned.insert(*it);
    }

    // Create object
    // Warning FP: comment next line for fast debug. Review Interaction construction: first Nslaw and Relation and THEN interaction
    // with nslaw and relation as arguments of the constructor list.
    //Interaction *interaction = new  Interaction(name,dsConcerned, nId,nbRel );    // interaction->display();

    // Push the interaction on global DS vectors
    //GLOB_SET_INTERACTION.insert(interaction);

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
  int nId = 0, dimDS = 0, nbRel = 0;

  try
  {
    InteractionsSet * InteractionSetPtr = GLOB_DATA.getInteractionsSetPtr();

    // Retrieve Interaction
    int nSize = InteractionSetPtr->size();
    if ((nIdInteraction < 0) || (nIdInteraction > nSize))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR failed");

    Interaction *interaction = InteractionSetPtr->getInteraction(nIdInteraction);

    // Compute sum of DS size
    dimDS = interaction->getSizeOfDS();
    nbRel = interaction->getSizeOfY();

    if ((dimDS <= 0) || (nbRel <= 0))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR  failed");

    // Vectors and Matrix creation
    SimpleMatrix mH(nbRel, dimDS);
    SimpleVector  vb(nbRel);

    // Vectors and Matrix initialisation with function parameters
    // Is there a solution less stupid ?
    //vb(0,0);
    for (int i = 0; i < nbRel; i++)
    {
      vb(i) = b[i];
      for (int j = 0; j < dimDS; j++)
      {
        mH(i, j) = H[i * nbRel + j];
      }
    }
    //    cout <<"H("<<dimDS<<")="<< mH <<endl;

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

/*
extern "C" int sicLagrangianR(int nIdInteraction, char *relationType, char *funcH, char *funcG)
{
  int nId=SIC_OK;

  try {
    InteractionsSet * InteractionSetPtr=GLOB_DATA.getInteractionsSetPtr();

    // Retrieve Interaction
    int nSize=InteractionSetPtr->size();
    if ( (nIdInteraction<0) || (nIdInteraction>nSize) )
      RuntimeException::selfThrow("siconos/C:: sicLagrangianR failed");

    Interaction *interaction = InteractionSetPtr->getInteraction(nIdInteraction);

    if ( interaction==NULL )
      RuntimeException::selfThrow("siconos/C:: sicLagrangianR failed");

    // Why a vector ? See Franck or Vincent...
    vector<string> vG;
    vG.reserve(1);
    vG.push_back(funcG);

    Relation * relation = new LagrangianR(relationType, funcH,vG);
    interaction->setRelationPtr(relation);
  }
  catch(SiconosException e)
    {cout << e.report() << endl;nId=SIC_ERROR;}
  catch(...)
    {cout << "Exception caught in sicLagrangianLinearR" << endl;nId=SIC_ERROR;}

  return nId;
}
*/

extern "C" int sicNewtonImpactNSL(int nIdInteraction, double e)
{
  int nId = SIC_OK;

  try
  {
    InteractionsSet * InteractionSetPtr = GLOB_DATA.getInteractionsSetPtr();

    // Retrieve Interaction
    int nSize = InteractionSetPtr->size();
    if ((nIdInteraction < 0) || (nIdInteraction > nSize))
      RuntimeException::selfThrow("siconos/C:: sicLagrangianLinearR failed");

    Interaction *interaction = InteractionSetPtr->getInteraction(nIdInteraction);

    NonSmoothLaw *nslaw = new NewtonImpactNSL(e);

    interaction->setNonSmoothLawPtr(nslaw);

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in ssicNewtonImpactNSL" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}


extern "C" int sicNonSmoothDynamicalSystem(int isBVP)
{
  int nId = SIC_OK;

  try
  {
    NonSmoothDynamicalSystem *prtNSDS = GLOB_DATA.getNonSmoothDynamicalSystemPtr();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();
    InteractionsSet * InteractionSetPtr = GLOB_DATA.getInteractionsSetPtr();

    if ((isBVP < 0) || (isBVP > 1))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad isBVP");
    // Create NSDS system connected with DS and Interactions
    // Warning FP: NSDS constructor change -> add set of DS and set of Interactions in arguments list.
    // => To be reviewed?
    prtNSDS = new NonSmoothDynamicalSystem(*SetDSPtr, *InteractionSetPtr, isBVP);

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
    Model * prtModel = GLOB_DATA.getModelPtr();
    NonSmoothDynamicalSystem *prtNSDS = GLOB_DATA.getNonSmoothDynamicalSystemPtr();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();
    InteractionsSet * InteractionSetPtr = GLOB_DATA.getInteractionsSetPtr();


    // Parameters verification
    if ((t0 < 0.0) || (T < 0.0))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad time parameters");

    if (SetDSPtr->size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to DS empty");


    if (InteractionSetPtr->size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to Interaction empty");
    if (prtNSDS == NULL)
      RuntimeException::selfThrow("siconos/C:: sicModel failed due to NSDS empty");

    prtNSDS->setDynamicalSystems(*SetDSPtr);
    prtNSDS->setInteractions(*InteractionSetPtr);

    // Create the model connected with NSDS
    prtModel = new Model(t0, T);
    prtModel->setNonSmoothDynamicalSystemPtr(prtNSDS);
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

extern "C" int sicSimulationTimeStepping(double h)
{
  int nId = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    Model * prtModel = GLOB_DATA.getModelPtr();
    TimeDiscretisation *prtTime = GLOB_DATA.getTimeDiscretisationPtr();

    prtSimul = new TimeStepping(prtTime);

    // Time discretisation
    prtTime = new  TimeDiscretisation(h, prtModel);

    prtSimul->initialize();

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSimulationTimeStepping" << endl;
    nId = SIC_ERROR;
  }

  return nId;

}


extern "C" int sicOneStepIntegratorMoreau(double *theta)
{
  int nId = SIC_OK, i;

  try
  {
    Simulation * prtSimul = GLOB_DATA.getSimulationPtr();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.getDynamicalSystemsSetPtr();

    int  dsNumber = SetDSPtr->size();

    if (dsNumber == 0)
      RuntimeException::selfThrow("siconos/C:: ssicSimulationTimeStepping failed due to DS empty");

    // One Step integrator s
    set<OneStepIntegrator *> vOSI;
    DSIterator it;
    i = 0;
    // \Warning Franck: corrections = consequences of Simulation/NSDS changes (vector<> => set<> )
    // Thus this part has to be reviewed -> especially the way theta values are sorted?
    for (it = SetDSPtr->begin(); it != SetDSPtr->end(); ++it)
    {
      vOSI.insert(new Moreau(*it, theta[i], prtSimul));
      i++;
    }

    prtSimul->setOneStepIntegrators(vOSI);
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
    //OneStepNSProblem * onsspb = new LCP(GLOB_SIMULATION,solverName,maxiter,tolerance);
    //    GLOB_SIMULATION->setOneStepNSProblemPtr(onsspb);

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
  //delete GLOB_MODEL;
  //delete GLOB_NSDS;
  //delete GLOB_SIMULATION;

  /*
  DSIterator itDS;
  for(itDS=GLOB_SET_DS->begin();itDS!=GLOB_SET_DS->end();++itDS)
    if( (*itDS) != NULL ) delete *itDS;

  InteractionsIterator it;
  for(it=GLOB_SET_INTERACTION.begin();it!=GLOB_SET_INTERACTION.end();++it)
    if ((*it) != NULL ) delete *it;
  */


  return SIC_OK;
}
