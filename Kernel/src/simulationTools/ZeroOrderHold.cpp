/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "ZeroOrderHold.hpp"
#include "EventDriven.hpp"
#include "Model.hpp"
#include "Lsodar.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "SubPluggedObject.hpp"
#include "FirstOrderType2R.hpp"
#include "MatrixIntegrator.hpp"

//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>


using namespace RELATION;
// --- constructor from a set of data ---
//ZeroOrderHold::ZeroOrderHold():
//  OneStepIntegrator(OSI::ZOH)
//{
//}

// --- xml constructor ---
ZeroOrderHold::ZeroOrderHold(SP::OneStepIntegratorXML osiXML, SP::DynamicalSystemsSet dsList):
  OneStepIntegrator(OSI::ZOH), _useGammaForRelation(false)
{
  RuntimeException::selfThrow("ZeroOrderHold::xml constructor - not yet implemented.");
}

// --- constructor from a minimum set of data ---
ZeroOrderHold::ZeroOrderHold(SP::DynamicalSystem ds):
  OneStepIntegrator(OSI::ZOH), _useGammaForRelation(false)
{
  OSIDynamicalSystems->insert(ds);
}


// Note: OSIDynamicalSystems must disappear
void ZeroOrderHold::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}


//void ZeroOrderHold::setPhi(const SiconosMatrix& newPhi, SP::DynamicalSystem ds)
//{
//  if (!ds)
//    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - ds == NULL.");
//
//  // Check if ds is in the OSI
//  if (!OSIDynamicalSystems->isIn(ds))
//    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - ds does not belong to this Integrator ...");
//
//  unsigned int dsN = ds->number();
//  // Check dimensions consistency
//  unsigned int line = newPhi.size(0);
//  unsigned int col  = newPhi.size(1);
//
//  if (line != col) // Check that newValue is square
//    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal,ds) - newVal is not square! ");
//
//  unsigned int sizePhi = ds->getDim(); // n for first order systems, ndof for lagrangian.
//  if (line != sizePhi) // check consistency between newValue and dynamical system size
//    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal,ds) - unconsistent dimension between newVal and dynamical system to be integrated ");
//
//  // Memory allocation for W, if required
//  if (!_PhiMap[dsN]) // allocate a new W if required
//  {
//    _PhiMap[dsN].reset(new SimpleMatrix(newPhi));
//  }
//  else  // or fill-in an existing one if dimensions are consistent.
//  {
//    if (line == _PhiMap[dsN]->size(0) && col == _PhiMap[dsN]->size(1))
//      *(_PhiMap[dsN]) = newPhi;
//    else
//      RuntimeException::selfThrow("ZeroOrderHold::setPhi: inconsistent dimensions with problem size for given input matrix W");
//  }
//}
//
//void ZeroOrderHold::setPhiPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds)
//{
//  unsigned int line = newPtr->size(0);
//  unsigned int col  = newPtr->size(1);
//  if (line != col) // Check that newPtr is square
//    RuntimeException::selfThrow("ZeroOrderHold::setPhiPtr(newVal, ds) - newVal is not square! ");
//
//  if (!ds)
//    RuntimeException::selfThrow("ZeroOrderHold::setPhiPtr(newVal, ds) - ds == NULL.");
//
//  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
//  if (line != sizeW) // check consistency between newValue and dynamical system size
//    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - unconsistent dimension between newVal and dynamical system to be integrated ");
//
//  _PhiMap[ds->number()] = newPtr;                  // link with new pointer
//}




void ZeroOrderHold::initialize()
{
  OneStepIntegrator::initialize();
  // Get initial time
  ConstDSIterator itDS;
  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  InteractionsGraph& IG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet0();
  DynamicalSystemsGraph::OEIterator oei, oeiend;
  Type::Siconos dsType;
  unsigned int indxIter;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    dsType = Type::value(**itDS);
    if ((dsType != Type::FirstOrderLinearDS) && (dsType != Type::FirstOrderLinearTIDS))
      RuntimeException::selfThrow("ZeroOrderHold::initialize - the DynamicalSystem does not have the right type");
    indxIter = 0;
    DynamicalSystemsGraph::AVIterator avi, aviend;
    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(*itDS);
    if (!DSG0.Ad.hasKey(dsgVD))
    {
      DSG0.Ad[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model()));
      if (DSG0.Ad.at(dsgVD)->isConst())
        DSG0.Ad.at(dsgVD)->integrate();
    }
    else
      RuntimeException::selfThrow("ZeroOrderHold::initialize - Ad MatrixIntegrator is already initialized for ds the DS");

    if ((static_cast<const FirstOrderLinearDS&>(**itDS)).b())
    {
      SP::SiconosMatrix E(new SimpleMatrix((*itDS)->getN(), (*itDS)->getN(), 0));
      E->eye();
      DSG0.AdInt.insert(dsgVD, SP::MatrixIntegrator(new MatrixIntegrator(**itDS, *simulationLink->model(), E)));
      if (DSG0.AdInt.at(dsgVD)->isConst())
        DSG0.AdInt.at(dsgVD)->integrate();
    }
   if (DSG0.B.hasKey(dsgVD))
   {
     DSG0.Bd[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), DSG0.B[dsgVD]));
     if (DSG0.Bd.at(dsgVD)->isConst())
       DSG0.Bd.at(dsgVD)->integrate();
   }
   if (DSG0.L.hasKey(dsgVD))
   {
     DSG0.Ld[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), DSG0.L[dsgVD]));
     if (DSG0.Ld.at(dsgVD)->isConst())
       DSG0.Ld.at(dsgVD)->integrate();
   }
   if (DSG0.pluginB.hasKey(dsgVD))
     DSG0.Bd[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), DSG0.pluginB[dsgVD], DSG0.u[dsgVD]->size()));
   if (DSG0.pluginL.hasKey(dsgVD))
     DSG0.Ld[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), DSG0.pluginL[dsgVD], DSG0.e[dsgVD]->size()));

    // Now we search for an Interaction dedicated to control
    for (std11::tie(avi, aviend) = DSG0.adjacent_vertices(dsgVD);
         avi != aviend; ++avi)
    {
      DynamicalSystemsGraph::EDescriptor ed1, ed2;
      std11::tie(ed1, ed2) = DSG0.edges(dsgVD, *avi);

      if (IG0.properties(IG0.descriptor(DSG0.bundle(ed1))).forControl)
      {
        Interaction& inter = *DSG0.bundle(ed1);
        Relation& rel = *inter.relation();
        if (rel.getType() != RELATION::FirstOrder)
          RuntimeException::selfThrow("ZeroOrderHold::initialize - the Integrator can only deal with FirstOrder Relation");
        FirstOrderR& relR = static_cast<FirstOrderR&>(rel);

        if (indxIter == 0)
        {
          indxIter++;
          if (!relR.isJacLgPlugged())
          {
            DSG0.Bd[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), relR.B()));
            if (DSG0.Bd.at(dsgVD)->isConst())
              DSG0.Bd.at(dsgVD)->integrate();
          }
          else
          {
            DSG0.Bd[dsgVD].reset(new MatrixIntegrator(**itDS, *simulationLink->model(), relR.getPluging(), inter.getSizeOfY()));
          }
        }
        else
        {
          //        RuntimeException::selfThrow("ZeroOrderHold::initialize - DS linked with more that one iteraction");
          DEBUG_PRINTF("number of iteraction attached to the process : %d\n", indxIter);
        }
      }
    }

    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, (*itDS)->getDim());

  }
}
//void ZeroOrderHold::initMatrixPhi(const double t, const DynamicalSystemsGraph::VDescriptor& dsgVD)
//{
//  // This function:
//  // - allocate memory for the matrices Phi matrix
//  // - insert this matrix into the DS properties
//
//  //  if (!ds)
//  //    RuntimeException::selfThrow("ZeroOrderHold::initMatrices(t,ds) - ds == NULL");
//
//
//  if (DSG0.Ad.hasKey(dsgVD))
//    RuntimeException::selfThrow("ZeroOrderHold::initMatrixPhi(t,ds) - initMatrixPhi(ds) is already in the map and has been initialized.");
//
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG0.bundle(dsgVD);
//  unsigned int sizeN = ds->getDim(); // n for first order systems
//
//  Type::Siconos dsType = Type::value(*ds);
//
//  // For first order linear time-invariant systems we can compute this once
//  // and for all
//  if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
//  {
//    DSG0.Ad[dsgVD].isConst = false;
//    DSG0.Ad[dsgVD].mat.reset(new SimpleMatrix(sizeN, sizeN));
//  }
//  else RuntimeException::selfThrow("ZeroOrderHold::initMatrixPhi - not yet implemented for Dynamical system type: " + dsType);
//
//}
//
//void ZeroOrderHold::initMatrices(const double t, SP::DynamicalSystem ds, const Interaction& inter)
//{
//  // This function:
//  // - allocate memory for the matrices Ad and Bd, Ld if needed
//
//  initMatrixPhi(t, ds);
//
//  unsigned int dsN = ds->number();
//  unsigned int sizeN = ds->getDim(); // n for first order systems
//  unsigned int sizeP = inter.getSizeOfY();         // p for first order systems
//
//  Type::Siconos dsType = Type::value(*ds);
//
//  // For first order linear time-invariant systems we can compute this once
//  // and for all
//  if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
//  {
//    _PsiMap[dsN].reset(new SimpleMatrix(sizeN, sizeP));
//    _isPsiConst[dsN] = false;
//
//  }
//  else RuntimeException::selfThrow("ZeroOrderHold::initMatrices - not yet implemented for Dynamical system type: " + dsType);
//
//}
//
//void ZeroOrderHold::initIntegrators(const DynamicalSystem& ds, const bool withInteraction)
//{
//  unsigned int dsN = ds.number();
//  double t0 = simulationLink->model()->t0();
//  double T = simulationLink->model()->finalT();
//  Type::Siconos dsType = Type::value(ds);
//
//  if (dsType == Type::FirstOrderLinearTIDS)
//  {
//    _DSPhiMap[dsN].reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
//    if (withInteraction)
//      _DSPsiMap[dsN].reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
//    _DSAdIntMap[dsN].reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
//  }
//  else if (dsType == Type::FirstOrderLinearDS)
//  {
//    const FirstOrderLinearDS& cfolds = static_cast<const FirstOrderLinearDS&>(ds);
//    _DSPhiMap[dsN].reset(new FirstOrderLinearDS(cfolds));
//    FirstOrderLinearDS& foldsPhi = static_cast<FirstOrderLinearDS&>(*_DSPhiMap[dsN]);
//    foldsPhi.zeroPlugin();
//    if (cfolds.getPluginA()->isPlugged())
//      foldsPhi.setPluginA(cfolds.getPluginA());
//     _DSAdIntMap[dsN].reset(new FirstOrderLinearDS(cfolds));
//    FirstOrderLinearDS& foldsAdInt = static_cast<FirstOrderLinearDS&>(*_DSAdIntMap[dsN]);
//    foldsAdInt.zeroPlugin();
//    if (cfolds.getPluginA()->isPlugged())
//      foldsAdInt.setPluginA(cfolds.getPluginA());
//   if (withInteraction)
//    {
//      _DSPsiMap[dsN].reset(new FirstOrderLinearDS(cfolds));
//      FirstOrderLinearDS& foldsPsi = static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]);
//      foldsPsi.zeroPlugin();
//      if (cfolds.getPluginA()->isPlugged())
//        foldsPsi.setPluginA(cfolds.getPluginA());
//    }
//  }
//
//  _modelPhiMap[dsN].reset(new Model(t0, T));
//  _modelPhiMap[dsN]->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPhiMap[dsN]);
//  _PhiOSIMap[dsN].reset(new Lsodar(_DSPhiMap[dsN]));
//  _TDPhiMap[dsN].reset(new TimeDiscretisation(*simulationLink->timeDiscretisation()));
//  _simulPhiMap[dsN].reset(new EventDriven(_TDPhiMap[dsN], 0));
//  _simulPhiMap[dsN]->insertIntegrator(_PhiOSIMap[dsN]);
//  _modelPhiMap[dsN]->initialize(_simulPhiMap[dsN]);
//
//  _modelAdIntMap[dsN].reset(new Model(t0, T));
//  _modelAdIntMap[dsN]->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSAdIntMap[dsN]);
//  _AdIntOSIMap[dsN].reset(new Lsodar(_DSAdIntMap[dsN]));
//  _TDAdIntMap[dsN].reset(new TimeDiscretisation(*simulationLink->timeDiscretisation()));
//  _simulAdIntMap[dsN].reset(new EventDriven(_TDAdIntMap[dsN], 0));
//  _simulAdIntMap[dsN]->insertIntegrator(_AdIntOSIMap[dsN]);
//  _modelAdIntMap[dsN]->initialize(_simulAdIntMap[dsN]);
//  if (withInteraction)
//  {
//    SP::SiconosVector tmpb(new SiconosVector(ds.getDim(), 0));
//    std11::static_pointer_cast<FirstOrderLinearDS>(_DSPsiMap[dsN])->setb(tmpb);
//    _modelPsiMap[dsN].reset(new Model(t0, T));
//    _modelPsiMap[dsN]->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPsiMap[dsN]);
//    _PsiOSIMap[dsN].reset(new Lsodar(_DSPsiMap[dsN]));
//    _TDPsiMap[dsN].reset(new TimeDiscretisation(*simulationLink->timeDiscretisation()));
//    _simulPsiMap[dsN].reset(new EventDriven(_TDPsiMap[dsN], 0));
//    _simulPsiMap[dsN]->insertIntegrator(_PsiOSIMap[dsN]);
//    _modelPsiMap[dsN]->initialize(_simulPsiMap[dsN]);
//  }
//}
//
//void ZeroOrderHold::computeAd(const DynamicalSystemsGraph::VDescriptor& dsgVD)
//{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG.bundle(dsgVD);
//  unsigned int dsN = ds.number();
//  unsigned int n = ds.getN();
//  SiconosVector* canonicalVector = new SiconosVector(n, 0);
//  EventDriven& sim = static_cast<EventDriven&>(*_simulPhiMap[dsN]);
//  SimpleMatrix& phi = DSG0.Ad[dsgVD].mat;
//  SiconosVector& x = *_DSPhiMap[dsN]->x();
//  //compute the matrix whose column are e^{A\Delta t}e_i
//  for (unsigned int i = 0; i < n; i++)
//  {
//    (*canonicalVector)(i) = 1;
//    x = *canonicalVector;
//    //Reset Lsodar
//    sim.setIstate(3);
//    sim.advanceToEvent();
//    phi.setCol(i, x);
//    (*canonicalVector)(i) = 0;
//  }
//  sim.processEvents();
//  delete canonicalVector;
//}
//
//void ZeroOrderHold::computeAdInt(const DynamicalSystemsGraph::VDescriptor& dsgVD)
//{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG.bundle(dsgVD);
//  unsigned int dsN = ds.number();
//  unsigned int n = ds.getN();
//  SiconosVector* canonicalVector = new SiconosVector(n, 0);
//  EventDriven& sim = static_cast<EventDriven&>(*_simulAdIntMap[dsN]);
//  SimpleMatrix& intPhi = *_AdIntMap[dsN];
//  SiconosVector& x = *_DSAdIntMap[dsN]->x();
//  SP::SiconosVector Bcol = static_cast<FirstOrderLinearDS&>(*_DSAdIntMap[dsN]).b();
//  if (!Bcol)
//  {
//    Bcol.reset(new SiconosVector(n, 0));
//    static_cast<FirstOrderLinearDS&>(*_DSAdIntMap[dsN]).setb(Bcol);
//  }
//
//  //compute the matrix whose column are \int e^{A\Delta \tau}d\tau e_i
//  for (unsigned int i = 0; i < n; i++)
//  {
//    (*canonicalVector)(i) = 1;
//    x.zero();
//    *Bcol = *canonicalVector;
//    //Reset Lsodar
//    sim.setIstate(3);
//    sim.advanceToEvent();
//    intPhi.setCol(i, x);
//    (*canonicalVector)(i) = 0;
//  }
//  sim.processEvents();
//  delete canonicalVector;
//}
//
//void ZeroOrderHold::computeBdTI(const DynamicalSystemsGraph::VDescriptor& dsgVD, const Relation& rel)
//{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG.bundle(dsgVD);
//  // Get relation and non smooth law types
//  RELATION::TYPES relationType = rel.getType();
//  RELATION::SUBTYPES relationSubType = rel.getSubType();
//  if ((relationType != FirstOrder) || (relationSubType != LinearTIR))
//    RuntimeException::selfThrow("ZeroOrderHold::computeBdTI - the associated Relation is not of type FirstOrderLinearTIR");
//  unsigned int dsN = ds.number();
//  EventDriven& sim = static_cast<EventDriven&>(*_simulPsiMap[dsN]);
//  SimpleMatrix& Bd = *DSG0.Bd[dsgVD];
//  SiconosVector& x = *_DSPsiMap[dsN]->x();
//  SiconosVector& Bcol = *static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]).b();
//  SiconosMatrix& B = *DSG0.B[dsgVD];
//  unsigned int p = B.size(1);
//  for (unsigned int i = 0; i < p; i++)
//  {
//    //compute the vector \int \Phi(t, \tau)\,\mathrm{d}\tau B_i
//    x.zero();
//    B.getCol(i, Bcol);
//    //Reset Lsodar
//    sim.setIstate(3);
//    //Compute
//    sim.processEvents();
//    sim.advanceToEvent();
//    Bd.setCol(i, x);
//  }
//}
//
//void ZeroOrderHold::computeLdTI(const DynamicalSystemsGraph::VDescriptor& dsgVD, const Relation& rel)
//{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG.bundle(dsgVD);
//  // Get relation and non smooth law types
//  unsigned int dsN = ds.number();
//  EventDriven& sim = static_cast<EventDriven&>(*_simulLdMap[dsN]);
//  SimpleMatrix& Ld = *DSG0.Ld[dsgVD];
//  SiconosVector& x = *_DSPsiMap[dsN]->x();
//  SiconosVector& Lcol = *static_cast<FirstOrderLinearDS&>(*_DSLdMap[dsN]).b();
//  SiconosMatrix& L = *DSG0.L[dsgVD];
//  unsigned int p = L.size(1);
//  for (unsigned int i = 0; i < p; i++)
//  {
//    //compute the vector \int \Phi(t, \tau)\,\mathrm{d}\tau B_i
//    x.zero();
//    L.getCol(i, Lcol);
//    //Reset Lsodar
//    sim.setIstate(3);
//    //Compute
//    sim.processEvents();
//    sim.advanceToEvent();
//    Ld.setCol(i, x);
//  }
//}
//
//void ZeroOrderHold::computeBd(const DynamicalSystemsGraph::VDescriptor& dsgVD, const Relation& rel)
//{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  DynamicalSystem& ds = DSG.bundle(dsgVD);
//  // Get relation and non smooth law types
//  RELATION::TYPES relationType = rel.getType();
//  RELATION::SUBTYPES relationSubType = rel.getSubType();
//  if ((relationType != FirstOrder) || (relationSubType != LinearTIR))
//    RuntimeException::selfThrow("ZeroOrderHold::computeBd - the associated Relation is not of type FirstOrderLinearTIR");
//  unsigned int dsN = ds.number();
//  EventDriven& sim = static_cast<EventDriven&>(*_simulPsiMap[dsN]);
//  SimpleMatrix& Bd = *DSG0.Bd[dsgVD];
//  FirstOrderLinearDS& foldsPsi = static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]);
//  SiconosVector& x = *foldsPsi.x();
//  SiconosVector& Bcol = *foldsPsi.b();
//  SiconosMatrix& B = *DSG0.B[dsgVD];
//  unsigned int p = B.size(1);
//  const FirstOrderLinearDS& folds = static_cast<const FirstOrderLinearDS&>(ds);
//  bool isAPlugged = folds.getPluginA()->isPlugged();
//  bool isBPlugged = rel.isJacLgPlugged();
//  if (isAPlugged)
//  {
//    RuntimeException::selfThrow("ZeroOrderHold::computeBd - Phi has to be constant for now");
//  }
//  SP::SubPluggedObject spo;
//  if (isBPlugged)
//  {
//    spo.reset(new SubPluggedObject(*rel.getPluginJacLg(), B.size(0), p));
//    foldsPsi.setPluginB(spo);
//  }
//  for (unsigned int i = 0; i < p; i++)
//  {
//    //compute the vector \int \Phi(t, \tau)\,\mathrm{d}\tau B_i
//    x.zero();
//    if (isBPlugged)
//      spo->setIndex(i);
//    else
//      B.getCol(i, Bcol);
//
//    //Reset Lsodar
//    sim.setIstate(3);
//    //Compute
//    sim.advanceToEvent();
//    Bd.setCol(i, x);
//  }
//  sim.processEvents();
//
//}

double ZeroOrderHold::computeResidu()
{

  // This function is used to compute the residu for each "Moreau-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it;
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it; // the considered dynamical system
    dsType = Type::value(ds); // Its type
    // 1 - First Order Linear Systems
    if (dsType == Type::FirstOrderLinearDS)
    {
      // No residu with ZOH ...
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      // No residu with ZOH ...
    }
    else
      RuntimeException::selfThrow("ZeroOrderHold::computeResidu - not yet implemented for Dynamical system type: " + dsType);
  }

  return maxResidu;
}

void ZeroOrderHold::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.

  DSIterator it; // Iterator through the set of DS.

  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  Type::Siconos dsType ; // Type of the current DS.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it; // the considered dynamical system
    dsType = Type::value(ds); // Its type

    DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(*it);
//    updateMatrices(dsDescr);
    if (dsType == Type::FirstOrderLinearTIDS || dsType == Type::FirstOrderLinearDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(ds);
      // Check whether we have to recompute things
      if (!DSG0.Ad.at(dsgVD)->isConst())
        DSG0.Ad.at(dsgVD)->integrate();
      if (d.b() && !DSG0.AdInt.at(dsgVD)->isConst())
        DSG0.AdInt.at(dsgVD)->integrate();
      if (DSG0.Bd.hasKey(dsgVD) && !DSG0.Bd.at(dsgVD)->isConst())
        DSG0.Bd.at(dsgVD)->integrate();
      if (DSG0.Ld.hasKey(dsgVD) && !DSG0.Ld.at(dsgVD)->isConst())
        DSG0.Ld.at(dsgVD)->integrate();

      SiconosVector& xfree = *d.workspace(DynamicalSystem::free);
      prod(DSG0.Ad.at(dsgVD)->mat(), *d.x(), xfree); // xfree = Ad*xold
      if (d.b())
      {
        assert(DSG0.AdInt.hasKey(dsgVD));
        prod(DSG0.AdInt.at(dsgVD)->mat(), *d.b(), xfree, false); // xfree += AdInt*b
      }
      // check whether we have a system with a control input
      if (DSG0.u.hasKey(dsgVD))
      {
        if (!DSG0.Bd.hasKey(dsgVD))
        {
          assert(DSG0.B.hasKey(dsgVD));
          DSG0.Bd[dsgVD].reset(new MatrixIntegrator(d, *simulationLink->model(), DSG0.B.at(dsgVD)));
          if (DSG0.Bd.at(dsgVD)->isConst())
            DSG0.Bd.at(dsgVD)->integrate();
        }
        prod(DSG0.Bd.at(dsgVD)->mat(), *DSG0.u.at(dsgVD), xfree, false); // xfree += Bd*u
      }
      // check whether the DynamicalSystem is an Observer
      if (DSG0.e.hasKey(dsgVD))
      {
        if (!DSG0.Ld.hasKey(dsgVD))
        {
          assert(DSG0.L.hasKey(dsgVD));
          DSG0.Ld[dsgVD].reset(new MatrixIntegrator(d, *simulationLink->model(), DSG0.L.at(dsgVD)));
          if (DSG0.Ld.at(dsgVD)->isConst())
            DSG0.Ld.at(dsgVD)->integrate();
        }
        prod(DSG0.Ld.at(dsgVD)->mat(), *DSG0.e.at(dsgVD), xfree, false); // xfree += -Ld*e
      }

    }
    else
      RuntimeException::selfThrow("ZeroOrderHold::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }

}

void ZeroOrderHold::prepareNewtonIteration(double time)
{
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    //    computeMatrices(time, **itDS);
  }
}


struct ZeroOrderHold::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) :
    _osnsp(p), _inter(inter) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->getNonSmoothLawSize();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    subscal(e, *_inter->y_k(_osnsp->levelMin()), *(_inter->yp()), subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    (*_inter->yp())(0) +=  e * (*_inter->y_k(_osnsp->levelMin()))(0);

  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};


void ZeroOrderHold::computeFreeOutput(SP::Interaction inter, OneStepNSProblem * osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
  */


  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->getRelationType();
  RELATION::SUBTYPES relationSubType = inter->getRelationSubType();

  unsigned int sizeY = inter->getNonSmoothLawSize();

  unsigned int relativePosition = 0;



  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;


  // All of these values should be stored in the node corrseponding to the UR when a Moreau scheme is used.
  SP::BlockVector Xq;
  Xq = inter->dataXq();
  SP::SiconosVector Yp = inter->yp();

  SP::BlockVector Xfree;
  if (relationType == FirstOrder)
  {
    Xfree = inter->data(FirstOrderR::free);
  }
  else if (relationType == NewtonEuler)
  {
    Xfree = inter->data(NewtonEulerR::free);
  }
  else if (relationType == Lagrangian)
  {
    Xfree = inter->data(LagrangianR::free);
  }
  // else if  (relationType == NewtonEuler)
  // {
  //   Xfree = inter->data(NewtonEulerR::free);
  // }
  assert(Xfree);

  SP::Relation rel = inter->relation();
  assert(inter);
  assert(rel);

  //  if (!IG0.properties(IG0.descriptor(inter)).forControl) // the integration is not for control
  {
    if (relationType == FirstOrder && relationSubType == Type2R)
    {
      SP::SiconosVector lambda = inter->lambda(0);
      SP::SiconosMatrix C = rel->C();
      SP::SiconosMatrix D = std11::static_pointer_cast<FirstOrderType2R>(rel)->D();
      assert(lambda);

      if (D)
      {
        coord[3] = D->size(1);
        coord[5] = D->size(1);
        subprod(*D, *lambda, *Yp, coord, true);

        *Yp *= -1.0;
      }
      if (C)
      {
        coord[3] = C->size(1);
        coord[5] = C->size(1);
        subprod(*C, *Xq, *Yp, coord, false);

      }

      if (_useGammaForRelation)
      {
        RuntimeException::selfThrow("ZeroOrderHold::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
      }
      SP::SiconosVector H_alpha = inter->Halpha();
      assert(H_alpha);
      *Yp += *H_alpha;
    }

    else
    {
      SP::SiconosMatrix C = rel->C();

      if (C)
      {

        assert(Xfree);
        assert(Xq);

        coord[3] = C->size(1);
        coord[5] = C->size(1);
        // creates a POINTER link between workX[ds] (xfree) and the
        // corresponding interactionBlock in each Interactionfor each ds of the
        // current Interaction.

        if (_useGammaForRelation)
        {
          subprod(*C, *Xq, *Yp, coord, true);
        }
        else
        {
          subprod(*C, *Xfree, *Yp, coord, true);
        }
      }

      if (relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
      {
        // In the first order linear case it may be required to add e + FZ to q.
        // q = HXfree + e + FZ
        SP::SiconosVector e;
        SP::SiconosMatrix F;
        if (relationSubType == LinearTIR)
        {
          e = std11::static_pointer_cast<FirstOrderLinearTIR>(rel)->e();
          F = std11::static_pointer_cast<FirstOrderLinearTIR>(rel)->F();
        }
        else
        {
          e = std11::static_pointer_cast<FirstOrderLinearR>(rel)->e();
          F = std11::static_pointer_cast<FirstOrderLinearR>(rel)->F();
        }

        if (e)
          *Yp += *e;

        if (F)
        {
          coord[3] = F->size(1);
          coord[5] = F->size(1);
          subprod(*F, *inter->dataZ(), *Yp, coord, false);
        }
      }

    }
  }


}
void ZeroOrderHold::integrate(double& tinit, double& tend, double& tout, int&)
{
  // This function should not be used
  RuntimeException::selfThrow("ZeroOrderHold::integrate - should not be used");
}

void ZeroOrderHold::updateState(const unsigned int level)
{
  bool useRCC = simulationLink->useRelativeConvergenceCriteron();
  if (useRCC)
    simulationLink->setRelativeConvergenceCriterionHeld(true);

  DSIterator it;
  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::OEIterator oei, oeiend;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it;
    Type::Siconos dsType = Type::value(ds);
    SP::Interaction interC;

    // 1 - First Order Linear Systems
    if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(ds);
      SiconosVector& x = *d.x();
      // 1 - First Order Linear Time Invariant Systems
      // \Phi is already computed
      x = *d.workspace(DynamicalSystem::free); // x = xfree = Phi*xold (+ Bd*u ) (+  Ld*e)
      if (level != LEVELMAX)
      {
        DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(*it);
       // we have to find the control interaction
       for (std11::tie(oei, oeiend) = DSG0.out_edges(dsgVD); oei != oeiend; ++oei)
       {
         if (DSG0.properties(*oei).forControl)
         {
           interC = DSG0.bundle(*oei);
           break;
         }
       }
        if (interC)
        {
          MatrixIntegrator& Bd = *DSG0.Bd[dsgVD];
          if (!Bd.isConst())
          {
            Bd.integrate();
          }
          prod(Bd.mat(), *interC->lambda(0), x, false); // x += Bd*\lambda
        }
      }
    }
    else
      RuntimeException::selfThrow("ZeroOrderHold::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


bool ZeroOrderHold::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF("ZeroOrderHold::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
#if __cplusplus >= 201103L
  assert(!::isnan(y));
#else
  assert(!std::isnan(y));
#endif
  if (y <= 0)
  {
    DEBUG_PRINT("ZeroOrderHold::addInteractionInIndexSet ACTIVATE.\n");
  }
  return (y <= 0);
}


bool ZeroOrderHold::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
  double gamma = .5;
  DEBUG_PRINTF("ZeroOrderHold::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
  y += gamma * h * yDot;
#if __cplusplus >= 201103L
  assert(!::isnan(y));
#else
  assert(!std::isnan(y));
#endif
  if (y > 0)
  {
    DEBUG_PRINT("ZeroOrderHold::removeInteractionInIndexSet DEACTIVATE.\n");
  }
  return (y > 0);
}

const SiconosMatrix& ZeroOrderHold::Ad(SP::DynamicalSystem ds)
{
  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  return DSG0.Ad.at(dsgVD)->mat();
}

const SiconosMatrix& ZeroOrderHold::Bd(SP::DynamicalSystem ds)
{
  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(ds);
  return DSG0.Bd.at(dsgVD)->mat();
}


void ZeroOrderHold::display()
{
  OneStepIntegrator::display();

  std::cout << "====== ZOH OSI display ======" <<std::endl;
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    std::cout << "--------------------------------" << std::endl;
//    cout << "--> Phi of dynamical system number " << itN << ": " << endl;
//    if (Ad(*it)) Ad(*it)->display();
//    else cout << "-> NULL" << endl;
//    cout << "--> Psi of dynamical system number " << itN << ": " << endl;
//    if (Bd(*it)) Bd(*it)->display();
//    else cout << "-> NULL" << endl;
  }
  std::cout << "================================" <<std::endl;
}

ZeroOrderHold* ZeroOrderHold::convert(OneStepIntegrator* osi)
{
  ZeroOrderHold* zeroOrderHold = dynamic_cast<ZeroOrderHold*>(osi);
  return zeroOrderHold;
}

void ZeroOrderHold::updateMatrices(SP::DynamicalSystem ds)
{
//  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
//  if (!DSG0.Ad[dsgVD]->isConst())
//    computeAd(dsgVD);
//  if (DSG0.Bd.hasKey(dsgVD) && !DSG->Bd[dsgVD]->isConst())
//    computeBd(ds);
}
