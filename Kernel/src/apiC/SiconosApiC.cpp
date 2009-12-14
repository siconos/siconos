y/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

#include "SiconosKernel.hpp"
#include "SiconosDataC.hpp"
#include "SiconosApiC.hpp"

using namespace std;

DataC GLOB_DATA;

extern "C" int sicLoadModel(char ModelXmlFile[])
{
  int ret = SIC_OK;

  try
  {
    Model *ptrModel = new Model(ModelXmlFile) ;
    GLOB_DATA.setModelPtr(ptrModel);
    // shortcut of object pointers
    GLOB_DATA.setSimulationPtr(ptrModel->simulation());
    GLOB_DATA.setEventsManagerPtr(ptrModel->simulation()->eventsManager());

    GLOB_DATA.setStatus(DATAC_MODEL);
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

    if (GLOB_DATA.getStatus() != DATAC_MODEL)
      RuntimeException::selfThrow("ApiC:: MODEL must be constructed before sicInitSimulation");

    // TimeStepping* s = static_cast<TimeStepping*>(GLOB_DATA.simulation());

    Simulation *s = GLOB_DATA.simulation();

    cout << "-->" << GLOB_DATA.simulation() << endl;
    //GLOB_DATA.simulation()->initialize();
    s->initialize();
    cout << "END" << endl;

    GLOB_DATA.setStatus(DATAC_INIT);

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
    Simulation * prtSimul = GLOB_DATA.simulation();
    *H = prtSimul->timeDiscretisation()->geth();
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
    Simulation * prtSimul = GLOB_DATA.simulation();
    *N = prtSimul->timeDiscretisation()->getNSteps();
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

    GLOB_DATA.eventsManager()->processEvents();

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


extern "C" int sicSTSaveInMemory()
{
  int ret = SIC_OK;
  try
  {
    Simulation * prtSimul = GLOB_DATA.simulation();
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
    GLOB_DATA.simulation()->advanceToEvent();
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
    Simulation * prtSimul = GLOB_DATA.simulation();
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
    Simulation * prtSimul = GLOB_DATA.simulation();
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

extern "C" int sicAdvanceToEvent()
{
  int ret = SIC_OK;
  try
  {

    cout << "sicAdvanceToEvent APIC" << endl;
    cout << "TYPE " << ((EventDriven*)GLOB_DATA.simulation())->getType() << endl;
    cout << "OK" << endl;
    ((EventDriven*) GLOB_DATA.simulation())->advanceToEvent();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicAdvanceToEvent" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicProcessEvents()
{
  int ret = SIC_OK;
  try
  {
    // Simulation * prtSimul=GLOB_DATA.simulation();
    //RPG OUPS !!!
    GLOB_DATA.eventsManager()->processEvents();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicProcessEvents" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicHasNextEvent(int *hasnextevent)
{
  int ret = SIC_OK;

  try
  {
    Simulation * prtSimul = GLOB_DATA.simulation();

    *hasnextevent = (int) prtSimul->hasNextEvent();
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicProcessEvents" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}

extern "C" int sicGetTypeEvent(char *type)
{
  int ret = SIC_OK;
  string typeSring;

  try
  {
    EventsManager *prtEventsManager = GLOB_DATA.eventsManager();
    strcpy(type, prtEventsManager->startingEvent()->getType().c_str());
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    ret = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicProcessEvents" << endl;
    ret = SIC_ERROR;
  }

  return ret;
}



extern "C" void sicDebug(int *ret)
{
  // GLOB_SIMULATION->update();
  cout << "--- Debug" << endl;

  // LagrangianDS* ball1 = static_cast<LagrangianDS*> (GLOB_MODEL->nonSmoothDynamicalSystem()->dynamicalSystem(0));

  //   vector<Interaction*> vIS= (GLOB_MODEL->nonSmoothDynamicalSystem()->getInteractions());

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

    Model * prtModel = GLOB_DATA.model();

    LagrangianDS* system = static_cast<LagrangianDS*>(prtModel->nonSmoothDynamicalSystem()->dynamicalSystem(indexDS));

    if (system)
    {
      size = (system->getQ()).size();
      if ((indexVector >= size) || (indexVector < 0))
        RuntimeException::selfThrow("siconos/C:: sicModelgetQ failed");
      else
      {
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


extern "C" int sicLagrangianLinearTIDS(int nDof, double *Q0, double *Vel0, double *Mass)
{
  int nId, i, j;

  try
  {

    DynamicalSystemsSet* SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // Id of LagrangianLinearTIDS
    nId = SetDSPtr->size();

    // TODO: parameters verification

    // Vectors and Matrix creation
    SimpleVector vQ0(nDof);
    SimpleVector vVel0(nDof);
    SimpleMatrix  mMass(nDof, nDof);


    // Vectors and Matrix initialisation with function parameters
    // Is there a solution less stupid ?
    for (i = 0; i < nDof; i++)
    {
      vQ0(i) = Q0[i];
      vVel0(i) = Vel0[i];
      for (j = 0; j < nDof; j++)
      {
        mMass(i, j) = Mass[j + i * nDof];
      }
    }

    // Push the LagrangianLinearTIDS on global DS vectors
    SetDSPtr->insert(new LagrangianLinearTIDS(nId, vQ0, vVel0, mMass));

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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

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
      DS->setComputeJacobianqFIntFunction(libname,fctJacFInt);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }

    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeMassFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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

extern "C" int sicSetComputeJacobianqNNLFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianqNNLFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianqNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)->setComputeJacobianNNLFunction(0, libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeJacobianqNNLFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int  sicSetComputeJacobianVelocityNNLFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C::  sicSetComputeJacobianVelocityNNLFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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

extern "C" int sicSetComputeJacobianqFIntFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianqFIntFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianqFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
    static_cast<LagrangianDS*>(DS)->setComputeJacobianFIntFunction(0, libname, func);
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    st = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicSetComputeJacobianqFIntFunction" << endl;
    st = SIC_ERROR;
  }

  return st;
}

extern "C" int sicSetComputeJacobianVelocityFIntFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetComputeJacobianVelocityFIntFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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

extern "C" int  sicSetFExt(int nIdDs, double *tFext)
{
  int st = SIC_OK, ndof;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C:: sicSetFExt failed");
      st = SIC_ERROR;
    }

    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);

    ndof = DS->getDim();

    SP::SiconosVector  Vfext = new SimpleVector(ndof);
    for (int index = 0; index < ndof; index++)
      (*Vfext)(index) = tFext[index];

    static_cast<LagrangianDS*>(DS)->setFExtPtr(Vfext);

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

extern "C" int  sicSetComputeFExtFunction(int nIdDs, char *libname, char *func)
{
  int st = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    // DS verification
    int  nDSMax = SetDSPtr->size();
    if ((nIdDs > nDSMax) || (nIdDs < 0))
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    if (!SetDSPtr->getPtr(nIdDs))
    {
      RuntimeException::selfThrow("siconos/C::  setComputeFExtFunction failed");
      st = SIC_ERROR;
    }
    DynamicalSystem *DS = SetDSPtr->getPtr(nIdDs);
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

extern "C" int sicInteraction(char *name, int nbDS, int *DS, int idLaw, int idRelation, int nSize)
{
  int nId;

  try
  {
    DynamicalSystemsSet * DSSetPtr = GLOB_DATA.dynamicalSystemsSet();
    NonSmoothLawSet * NSLawSetPtr = GLOB_DATA.nonSmoothLawSet();
    RelationsSet *RelationsSetPtr = GLOB_DATA.relationsSet();
    InteractionsSet * InteractionsSetPtr = GLOB_DATA.interactionsSet();

    // Id of Interaction
    nId = InteractionsSetPtr->size();

    // TODO: parameters verification
    int  nDSMax = DSSetPtr->size();
    if ((nbDS < 0) || (nDSMax < nbDS))
      RuntimeException::selfThrow("siconos/C:: sicInteraction failed");

    int  nLawMax = NSLawSetPtr->size();
    if ((idLaw < 0) || (idLaw > nLawMax))
      RuntimeException::selfThrow("siconos/C:: sicInteraction failed");

    int  nRelMax = RelationsSetPtr->size();
    if ((idRelation < 0) || (idRelation > nRelMax))
      RuntimeException::selfThrow("siconos/C:: sicInteraction failed");

    // Compute sum of DS size and add DS into the set
    DSIterator it;
    DynamicalSystemsSet dsConcerned;
    for (it = DSSetPtr->begin(); it != DSSetPtr->end(); ++it)
    {
      //---- TODO ---- get ???
      //dimDS+=(*it)->getDim();
      dsConcerned.insert(*it);
    }

    // Interaction construction
    Relation *rel     = (*RelationsSetPtr)[idRelation];
    NonSmoothLaw *law = (*NSLawSetPtr)[idLaw];
    InteractionsSetPtr->insert(new Interaction(name, dsConcerned, nId, nSize, law, rel));

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



extern "C" int sicLagrangianLinearTIR(int nDof, int nRel, double *H, double *b)
{
  int nId = 0;

  try
  {
    RelationsSet * RelationsSetPtr = GLOB_DATA.relationsSet();

    // Id of sicLagrangianLinearTIR
    nId = RelationsSetPtr->size();

    SimpleMatrix mH(nRel, nDof);

    for (int i = 0; i < nRel; i++)
    {
      for (int j = 0; j < nDof; j++)
      {
        mH(i, j) = H[i * nRel + j];
      }
    }

    SimpleVector vB(nRel);

    for (int i = 0; i < nRel; i++)
    {
      vB(i) = b[i];
    }

    RelationsSetPtr->push_back(new LagrangianLinearTIR(mH, vB));
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicLagrangianLinearTIR" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}

/*
extern "C" int sicLagrangianR(int nIdInteraction, char *relationType, char *funcH, char *funcG)
{
  int nId=SIC_OK;

  try {
    InteractionsSet * InteractionSetPtr=GLOB_DATA.interactionsSet();

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
    {cout << "Exception caught in sicLagrangianLinearTIR" << endl;nId=SIC_ERROR;}

  return nId;
}
*/

extern "C" int sicNewtonImpactNSL(double e)
{
  int nId = 0;

  try
  {
    NonSmoothLawSet * NSLawSetPtr = GLOB_DATA.nonSmoothLawSet();

    // Id of NewtonImpactNSL
    nId = NSLawSetPtr->size();
    // Push the NewtonImpactNSL
    NSLawSetPtr->push_back(new NewtonImpactNSL(e));
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicNewtonImpactNSL" << endl;
    nId = SIC_ERROR;
  }

  return nId;
}


extern "C" int sicNonSmoothDynamicalSystem(int isBVP)
{
  int nId = SIC_OK;

  try
  {
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();
    InteractionsSet * InteractionSetPtr = GLOB_DATA.interactionsSet();

    if ((isBVP < 0) || (isBVP > 1))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad isBVP");

    GLOB_DATA.setNonSmoothDynamicalSystemPtr(new NonSmoothDynamicalSystem(*SetDSPtr, *InteractionSetPtr, isBVP));

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
    NonSmoothDynamicalSystem *prtNSDS = GLOB_DATA.nonSmoothDynamicalSystem();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();
    InteractionsSet * InteractionSetPtr = GLOB_DATA.interactionsSet();


    // Parameters verification
    if ((t0 < 0.0) || (T < 0.0))
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to bad time parameters");

    if (SetDSPtr->size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to DS empty");


    if (InteractionSetPtr->size() == 0)
      RuntimeException::selfThrow("siconos/C:: sicNSDSModel failed due to Interaction empty");
    if (!prtNSDS)
      RuntimeException::selfThrow("siconos/C:: sicModel failed due to NSDS empty");

    //prtNSDS->setDynamicalSystems(*SetDSPtr);
    //prtNSDS->setInteractions(*InteractionSetPtr);

    // Create the model connected with NSDS
    GLOB_DATA.setModelPtr(new Model(t0, T));
    GLOB_DATA.model()->setNonSmoothDynamicalSystemPtr(prtNSDS);
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

extern "C" int sicTimeDiscretisation(double h)
{
  int nId;

  try
  {
    Model * prtModel = GLOB_DATA.model();

    TimesSet * TimesSetPtr = GLOB_DATA.timesSet();

    // Id of  Time discretisation
    nId = TimesSetPtr->size();
    // Push the Time discretisation
    TimesSetPtr->push_back(new  TimeDiscretisation(h, prtModel));

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
    nId = SIC_ERROR;
  }
  catch (...)
  {
    cout << "Exception caught in sicTimeDiscretisation" << endl;
    nId = SIC_ERROR;
  }

  return nId;

}

extern "C" int sicSimulationTimeStepping(int idTime)
{
  int nId = SIC_OK;
  int nSize;

  try
  {
    // Time discretisation
    TimesSet * timediscret = GLOB_DATA.timesSet();
    nSize = timediscret->size();
    if ((idTime < 0) || (idTime > nSize))
      RuntimeException::selfThrow("siconos/C:: sicSimulationTimeStepping failed due to bad time index");

    GLOB_DATA.setSimulationPtr(new TimeStepping((*timediscret)[idTime]));

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
    Simulation * prtSimul = GLOB_DATA.simulation();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    int  dsNumber = SetDSPtr->size();

    if (dsNumber == 0)
      RuntimeException::selfThrow("siconos/C:: sicOneStepIntegratorMoreau failed due to DS empty");

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

/* TBV by Franck */
extern "C" int sicOneStepIntegratorLsodar()
{
  int nId = SIC_OK, i;

  try
  {
    Simulation * prtSimul = GLOB_DATA.simulation();
    DynamicalSystemsSet * SetDSPtr = GLOB_DATA.dynamicalSystemsSet();

    int  dsNumber = SetDSPtr->size();

    if (dsNumber == 0)
      RuntimeException::selfThrow("siconos/C:: sicOneStepIntegratorMoreau failed due to DS empty");

    // One Step integrator s
    set<OneStepIntegrator *> vOSI;
    DSIterator it;
    i = 0;
    // \Warning Franck: corrections = consequences of Simulation/NSDS changes (vector<> => set<> )
    // Thus this part has to be reviewed -> especially the way theta values are sorted?
    for (it = SetDSPtr->begin(); it != SetDSPtr->end(); ++it)
    {
      vOSI.insert(new Lsodar(*it, prtSimul));
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
    cout << "Exception caught in sicOneStepIntegratorLsodar" << endl;
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

    // TimeStepping *s = (TimeStepping *)GLOB_DATA.simulation();
    // SP::OneStepNSProblem onsspb = new LCP(GLOB_DATA.simulation(),"LCP",solverName,maxiter,tolerance,normetyp,searchdir,rho);
    // PB
    // SP::OneStepNSProblem onsspb = new LCP(s,"LCP",solverName,maxiter,tolerance);
    // FP: unused variables => comment.

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

  // shortcut of object pointers
  GLOB_DATA.setEventsManagerPtr(GLOB_DATA.simulation()->eventsManager());

  GLOB_DATA.setStatus(DATAC_MODEL);

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
