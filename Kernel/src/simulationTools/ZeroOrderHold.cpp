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

//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace std;
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

// --- constructor from a set of data ---
ZeroOrderHold::ZeroOrderHold(DynamicalSystemsSet& allDS):
  OneStepIntegrator(OSI::ZOH, allDS)
{
}

// --- constructor from a minimum set of data ---
ZeroOrderHold::ZeroOrderHold(SP::DynamicalSystem ds):
  OneStepIntegrator(OSI::ZOH)
{
  OSIDynamicalSystems->insert(ds);
}


// Note: OSIDynamicalSystems must disappear
void ZeroOrderHold::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}


void ZeroOrderHold::setPhi(const SiconosMatrix& newPhi, SP::DynamicalSystem ds)
{
  if (!ds)
    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - ds == NULL.");

  // Check if ds is in the OSI
  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - ds does not belong to this Integrator ...");

  unsigned int dsN = ds->number();
  // Check dimensions consistency
  unsigned int line = newPhi.size(0);
  unsigned int col  = newPhi.size(1);

  if (line != col) // Check that newValue is square
    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal,ds) - newVal is not square! ");

  unsigned int sizePhi = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizePhi) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal,ds) - unconsistent dimension between newVal and dynamical system to be integrated ");

  // Memory allocation for W, if required
  if (!_PhiMap[dsN]) // allocate a new W if required
  {
    _PhiMap[dsN].reset(new SimpleMatrix(newPhi));
  }
  else  // or fill-in an existing one if dimensions are consistent.
  {
    if (line == _PhiMap[dsN]->size(0) && col == _PhiMap[dsN]->size(1))
      *(_PhiMap[dsN]) = newPhi;
    else
      RuntimeException::selfThrow("ZeroOrderHold::setPhi: inconsistent dimensions with problem size for given input matrix W");
  }
}

void ZeroOrderHold::setPhiPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds)
{
  unsigned int line = newPtr->size(0);
  unsigned int col  = newPtr->size(1);
  if (line != col) // Check that newPtr is square
    RuntimeException::selfThrow("ZeroOrderHold::setPhiPtr(newVal, ds) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("ZeroOrderHold::setPhiPtr(newVal, ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("ZeroOrderHold::setPhi(newVal, ds) - unconsistent dimension between newVal and dynamical system to be integrated ");

  _PhiMap[ds->number()] = newPtr;                  // link with new pointer
}




void ZeroOrderHold::initialize()
{
  OneStepIntegrator::initialize();
  // Get initial time
  double t0 = simulationLink->model()->t0();
  // Initialize ans compute Phi and Bd for all ds
  ConstDSIterator itDS;
  DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
  InteractionsGraph& IG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->indexSet0();
  DynamicalSystemsGraph::OEIterator oei, oeiend;
  Type::Siconos dsType;
  unsigned int indxIter;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    dsType = Type::value(**itDS);
    indxIter = 0;
    DynamicalSystemsGraph::AVIterator avi, aviend;
    for (boost::tie(avi, aviend) = DSG0.adjacent_vertices(DSG0.descriptor(*itDS));
         avi != aviend; ++avi)
    {
      DynamicalSystemsGraph::EDescriptor ed1, ed2;
      boost::tie(ed1, ed2) = DSG0.edges(DSG0.descriptor(*itDS), *avi);

      if (IG0.properties(IG0.descriptor(DSG0.bundle(ed1))).forControl) // the integration is for control
      {
        Interaction& inter = *DSG0.bundle(ed1);

        Relation& rel = *inter.relation();

        if (indxIter == 0)
        {
          indxIter++;
          // Matrices and integrator initialization
          initMatrices(t0, *itDS, inter);
          initIntegrators(**itDS, true);
          if (dsType == Type::FirstOrderLinearTIDS)
          {
            FirstOrderLinearTIDS& d = *boost::static_pointer_cast<FirstOrderLinearTIDS>(*itDS);
            if (_constH)
            {
              computePhi(d);
              computePsiTI(d, rel);
            }
          }
          else if (dsType == Type::FirstOrderLinearDS)
          {
            FirstOrderLinearDS& d = *boost::static_pointer_cast<FirstOrderLinearDS>(*itDS);
            if (_constH)
            {
              computePhi(d);
              if (!(rel.isJacLgPlugged() || d.getPluginA()->isPlugged()))
              {
                // in fact everything is constant ...
                computePsiTI(d, rel);
              }
              else
                computePsi(d, rel);
            }
          }
        }
        else
          //        RuntimeException::selfThrow("ZeroOrderHold::initialize - DS linked with more that one iteraction");
          cout << "number of iteraction attached to the process : " << indxIter << endl;
      }
    }
    if (indxIter == 0)
    {
      initMatrixPhi(t0, *itDS);
      initIntegrators(**itDS, false);
      if (dsType == Type::FirstOrderLinearTIDS)
      {
        FirstOrderLinearTIDS& d = *boost::static_pointer_cast<FirstOrderLinearTIDS>(*itDS);
        if (_constH)
        {
          computePhi(d);
        }
      }
    }

    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, _PhiMap[(*itDS)->number()]->size(0));

  }
}
void ZeroOrderHold::initMatrixPhi(const double t, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for the matrices Phi matrix
  // - insert this matrix into WMap with ds as a key

  //  if (!ds)
  //    RuntimeException::selfThrow("ZeroOrderHold::initMatrices(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("ZeroOrderHold::initMatrixPhi(t,ds) - ds does not belong to the OSI.");

  unsigned int dsN = ds->number();

  if (_PhiMap.find(dsN) != _PhiMap.end())
    RuntimeException::selfThrow("ZeroOrderHold::initMatrixPhi(t,ds) - initMatrixPhi(ds) is already in the map and has been initialized.");

  unsigned int sizeN = ds->getDim(); // n for first order systems

  _constH = simulationLink->timeDiscretisation()->getTDCase() == 2 ? true : false;
  Type::Siconos dsType = Type::value(*ds);

  // For first order linear time-invariant systems we can compute this once
  // and for all
  if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    _PhiMap[dsN].reset(new SimpleMatrix(sizeN, sizeN));
    _xNext[dsN].reset(new SimpleVector(sizeN));
  }
  else RuntimeException::selfThrow("ZeroOrderHold::initMatrixPhi - not yet implemented for Dynamical system type: " + dsType);

}
void ZeroOrderHold::initMatrices(const double t, SP::DynamicalSystem ds, const Interaction& inter)
{
  // This function:
  // - allocate memory for the matrices Phi and Psi matrices
  // - insert this matrix into WMap with ds as a key

  //  if (!ds)
  //    RuntimeException::selfThrow("ZeroOrderHold::initMatrices(t,ds) - ds == NULL");

  initMatrixPhi(t, ds);

  unsigned int dsN = ds->number();
  unsigned int sizeN = ds->getDim(); // n for first order systems
  unsigned int sizeP = inter.getSizeOfY();         // p for first order systems

  _constH = simulationLink->timeDiscretisation()->getTDCase() == 2 ? true : false;
  Type::Siconos dsType = Type::value(*ds);

  // For first order linear time-invariant systems we can compute this once
  // and for all
  if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    _PsiMap[dsN].reset(new SimpleMatrix(sizeN, sizeP));

  }
  else RuntimeException::selfThrow("ZeroOrderHold::initMatrices - not yet implemented for Dynamical system type: " + dsType);

}

void ZeroOrderHold::initIntegrators(const DynamicalSystem& ds, const bool withInteraction)
{
  unsigned int dsN = ds.number();
  unsigned int t0 = simulationLink->model()->t0();
  unsigned int T = simulationLink->model()->finalT();
  Type::Siconos dsType = Type::value(ds);

  if (dsType == Type::FirstOrderLinearTIDS)
  {
    _DSPhiMap[dsN].reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
    if (withInteraction)
      _DSPsiMap[dsN].reset(new FirstOrderLinearTIDS(static_cast<const FirstOrderLinearTIDS&>(ds)));
    //    _DSPsiMap[dsN]->setx0(
  }
  else if (dsType == Type::FirstOrderLinearDS)
  {
    const FirstOrderLinearDS& cfolds = static_cast<const FirstOrderLinearDS&>(ds);
    _DSPhiMap[dsN].reset(new FirstOrderLinearDS(cfolds));
    FirstOrderLinearDS& foldsPhi = static_cast<FirstOrderLinearDS&>(*_DSPhiMap[dsN]);
    foldsPhi.zeroPlugin();
    if (cfolds.getPluginA()->isPlugged())
      foldsPhi.setPluginA(cfolds.getPluginA());
    if (withInteraction)
    {
      _DSPsiMap[dsN].reset(new FirstOrderLinearDS(cfolds));
      FirstOrderLinearDS& foldsPsi = static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]);
      foldsPsi.zeroPlugin();
      if (cfolds.getPluginA()->isPlugged())
        foldsPsi.setPluginA(cfolds.getPluginA());
    }
  }

  _modelPhiMap[dsN].reset(new Model(t0, T));
  _modelPhiMap[dsN]->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPhiMap[dsN]);
  _PhiOSIMap[dsN].reset(new Lsodar(_DSPhiMap[dsN]));
  _TDPhiMap[dsN].reset(new TimeDiscretisation(*simulationLink->timeDiscretisation()));
  _simulPhiMap[dsN].reset(new EventDriven(_TDPhiMap[dsN], 0));
  _simulPhiMap[dsN]->insertIntegrator(_PhiOSIMap[dsN]);
  _modelPhiMap[dsN]->initialize(_simulPhiMap[dsN]);

  if (withInteraction)
  {
    SP::SiconosVector tmpb(new SimpleVector(ds.getDim(), 0));
    static_pointer_cast<FirstOrderLinearDS>(_DSPsiMap[dsN])->setb(tmpb);
    _modelPsiMap[dsN].reset(new Model(t0, T));
    _modelPsiMap[dsN]->nonSmoothDynamicalSystem()->insertDynamicalSystem(_DSPsiMap[dsN]);
    _PsiOSIMap[dsN].reset(new Lsodar(_DSPsiMap[dsN]));
    _TDPsiMap[dsN].reset(new TimeDiscretisation(*simulationLink->timeDiscretisation()));
    _simulPsiMap[dsN].reset(new EventDriven(_TDPsiMap[dsN], 0));
    _simulPsiMap[dsN]->insertIntegrator(_PsiOSIMap[dsN]);
    _modelPsiMap[dsN]->initialize(_simulPsiMap[dsN]);
  }
}

void ZeroOrderHold::computeNextX(const DynamicalSystem& ds)
{
  unsigned int dsN = ds.number();
  EventDriven& sim = static_cast<EventDriven&>(*_simulPhiMap[dsN]);
  SiconosVector& x = *_DSPhiMap[dsN]->x();
  x = ds.getx();
  sim.setIstate(3);
  sim.processEvents();
  sim.advanceToEvent();
  *_xNext[dsN] = x;
}

void ZeroOrderHold::computePhi(const DynamicalSystem& ds)
{
  unsigned int dsN = ds.number();
  unsigned int n = ds.getN();
  SimpleVector* canonicalVector = new SimpleVector(n, 0);
  EventDriven& sim = static_cast<EventDriven&>(*_simulPhiMap[dsN]);
  SimpleMatrix& phi = *_PhiMap[dsN];
  SiconosVector& x = *_DSPhiMap[dsN]->x();
  //compute the matrix whose column are e^{A\Delta t}e_i
  for (unsigned int i = 0; i < n; i++)
  {
    (*canonicalVector)(i) = 1;
    x = *canonicalVector;
    //Reset Lsodar
    sim.setIstate(3);
    sim.processEvents();
    sim.advanceToEvent();
    phi.setCol(i, x);
    (*canonicalVector)(i) = 0;
  }
  delete canonicalVector;
}

void ZeroOrderHold::computePsiTI(const DynamicalSystem& ds, const Relation& rel)
{
  // Get relation and non smooth law types
  RELATION::TYPES relationType = rel.getType();
  RELATION::SUBTYPES relationSubType = rel.getSubType();
  if ((relationType != FirstOrder) || (relationSubType != LinearTIR))
    RuntimeException::selfThrow("ZeroOrderHold::computePsiTI - the associated Relation is not of type FirstOrderLinearTIR");
  unsigned int dsN = ds.number();
  EventDriven& sim = static_cast<EventDriven&>(*_simulPsiMap[dsN]);
  SimpleMatrix& psi = *_PsiMap[dsN];
  SiconosVector& x = *_DSPsiMap[dsN]->x();
  SiconosVector& Bcol = *static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]).b();
  SiconosMatrix& B = *static_cast<const FirstOrderLinearTIR&>(rel).B();
  unsigned int p = B.size(1);
  for (unsigned int i = 0; i < p; i++)
  {
    //compute the vector \int \Phi(t, \tau)\,\mathrm{d}\tau B_i
    x.zero();
    B.getCol(i, Bcol);
    //Reset Lsodar
    sim.setIstate(3);
    //Compute
    sim.processEvents();
    sim.advanceToEvent();
    psi.setCol(i, x);
  }
}

// XXX untested
void ZeroOrderHold::computePsi(const DynamicalSystem& ds, const Relation& rel)
{
  // Get relation and non smooth law types
  RELATION::TYPES relationType = rel.getType();
  RELATION::SUBTYPES relationSubType = rel.getSubType();
  if ((relationType != FirstOrder) || (relationSubType != LinearTIR))
    RuntimeException::selfThrow("ZeroOrderHold::computePsiTI - the associated Relation is not of type FirstOrderLinearTIR");
  unsigned int dsN = ds.number();
  EventDriven& sim = static_cast<EventDriven&>(*_simulPsiMap[dsN]);
  SimpleMatrix& psi = *_PsiMap[dsN];
  SimpleMatrix& phi = *_PhiMap[dsN];
  FirstOrderLinearDS& foldsPsi = static_cast<FirstOrderLinearDS&>(*_DSPsiMap[dsN]);
  SiconosVector& x = *foldsPsi.x();
  SiconosVector& Bcol = *foldsPsi.b();
  SiconosMatrix& B = *static_cast<const FirstOrderLinearTIR&>(rel).B();
  unsigned int p = B.size(1);
  const FirstOrderLinearDS& folds = static_cast<const FirstOrderLinearDS&>(ds);
  bool isAPlugged = folds.getPluginA()->isPlugged();
  bool isBPlugged = rel.isJacLgPlugged();
  if (isAPlugged)
  {
    RuntimeException::selfThrow("ZeroOrderHold::computePsiTI - Phi has to be constant for now");
  }
  else
    foldsPsi.setA(phi);
  SP::SubPluggedObject spo;
  if (isBPlugged)
  {
    spo.reset(new SubPluggedObject(*rel.getPluginJacLg(), B.size(0), p));
    foldsPsi.setPluginB(spo);
  }
  for (unsigned int i = 0; i < p; i++)
  {
    //compute the vector \int \Phi(t, \tau)\,\mathrm{d}\tau B_i
    x.zero();
    if (isBPlugged)
      spo->setIndex(i);
    else
      B.getCol(i, Bcol);

    //Reset Lsodar
    sim.setIstate(3);
    //Compute
    sim.processEvents();
    sim.advanceToEvent();
    psi.setCol(i, x);
  }

}

void ZeroOrderHold::computeMatrices(const double t, const DynamicalSystem& ds)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, WMap[ds] is supposed to exist and not to be null
  // Memory allocation has been done during initW.

  //assert(ds && "ZeroOrderHold::computeMatrices(t,ds) - ds == NULL");

  assert((_PhiMap.find(ds.number()) != _PhiMap.end()) &&
         "ZeroOrderHold::computeMatrices(t,ds) - Phi(ds) does not exists. Maybe you forget to initialize the osi?");

  Type::Siconos dsType = Type::value(ds);


  // 1 - First order linear TI systems
  if (dsType == Type::FirstOrderLinearTIDS)
  {
    //Nothing to be done here
  }
  // 2 - First order linear systems
  else if (dsType == Type::FirstOrderLinearDS)
  {
  }
  else RuntimeException::selfThrow("ZeroOrderHold::computeMatrices - not yet implemented for Dynamical system type :" + dsType);
}




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
  double normResidu = maxResidu;

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it; // the considered dynamical system
    dsType = Type::value(ds); // Its type
    SiconosVector& residuFree = *ds.residuFree();
    // 1 - First Order Linear Systems
    if (dsType == Type::FirstOrderLinearDS)
    {
      // No residu with ZOH ...
      residuFree.zero();
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      // No residu with ZOH ...
      residuFree.zero();
    }
    else
      RuntimeException::selfThrow("ZeroOrderHold::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void ZeroOrderHold::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  // Operators computed at told have index i, and (i+1) at t.

  DSIterator it; // Iterator through the set of DS.

  Type::Siconos dsType ; // Type of the current DS.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it; // the considered dynamical system
    dsType = Type::value(ds); // Its type
    unsigned int dsN = ds.number();
    // 1 - First Order Linear Time Invariant Systems
    // The transition matrix is constant
    // if the timestep stays the same the computation is really easy
    if (dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearTIDS& d = static_cast<FirstOrderLinearTIDS&>(ds);
      SiconosVector& xfree = *d.workFree();//workX[d];
      SiconosMatrix& Phi = *_PhiMap[dsN]; // Phi is constant
      prod(Phi, *d.x(), xfree); // done
      *_xNext[dsN] = xfree;
    }
    else if (dsType == Type::FirstOrderLinearDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(ds);
      SiconosVector& xfree = *(d.workFree());

      // Update _xNext
      //computeNextX(ds);
      computePhi(d);
      SiconosMatrix& Phi = *_PhiMap[dsN];
      prod(Phi, *d.x(), xfree);
      *_xNext[dsN] = xfree;
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


  OneStepNSProblems&  allOSNS = *simulationLink->oneStepNSProblems();

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
  SP::SiconosVector Xq = inter->workXq();
  SP::SiconosVector Yp = inter->yp();

  SP::SiconosVector Xfree = inter->workFree();

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
      SP::SiconosMatrix D = rel->D();
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
      SP::SiconosVector H_alpha = rel->Halpha();
      assert(H_alpha);
      *Yp += *H_alpha;
    }

    else if (relationType == NewtonEuler)
    {
      SP::SiconosMatrix CT =  static_pointer_cast<NewtonEulerR>(rel)->jachqT();

      if (CT)
      {

        assert(Xfree);
        if (Xfree->size() == 0)
        {
          // creates a POINTER link between workX[ds] (xfree) and the
          // corresponding unitaryBlock in each UR for each ds of the
          // current UR.
          ConstDSIterator itDS;
          for (itDS = inter->dynamicalSystemsBegin();
               itDS != inter->dynamicalSystemsEnd();
               ++itDS)
          {
            //osi = osiMap[*itDS];
            inter->insertInWorkFree((*itDS)->workFree()); // osi->getWorkX(*itDS));
          }
        }
        assert(Yp);

        coord[3] = CT->size(1);
        coord[5] = CT->size(1);
        // printf("LinearOSNS: computing q: CT\n");
        // CT->display();
        // printf("LinearOSNS: computing q: Xfree and _q\n");
        // Xfree->display();
        subprod(*CT, *Xfree, *Yp, coord, true);
        //        _q->display();
      }

    }
    else
    {
      SP::SiconosMatrix C = rel->C();

      if (C)
      {

        assert(Xfree);
        if (Xfree->size() == 0)
        {
          // creates a POINTER link between workX[ds] (xfree) and the
          // corresponding unitaryBlock in each UR for each ds of the
          // current UR.
          ConstDSIterator itDS;
          for (itDS = inter->dynamicalSystemsBegin();
               itDS != inter->dynamicalSystemsEnd();
               ++itDS)
          {
            //osi = osiMap[*itDS];
            inter->insertInWorkFree((*itDS)->workFree()); // osi->getWorkX(*itDS));
          }
        }
        assert(Yp);
        assert(Xq);

        coord[3] = C->size(1);
        coord[5] = C->size(1);
        if (_useGammaForRelation)
        {
          subprod(*C, *Xq, *Yp, coord, true);
        }
        else
        {
          subprod(*C, *Xfree, *Yp, coord, true);
        }
      }

      if (relationType == Lagrangian)
      {
        SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
        ID->eye();

        Index xcoord(8);
        xcoord[0] = 0;
        xcoord[1] = sizeY;
        xcoord[2] = 0;
        xcoord[3] = sizeY;
        xcoord[4] = 0;
        xcoord[5] = sizeY;
        xcoord[6] = 0;
        xcoord[7] = sizeY;

        // For the relation of type LagrangianRheonomousR
        if (relationSubType == RheonomousR)
        {
          if (((allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
          {
            static_pointer_cast<LagrangianRheonomousR>(rel)->computehDot(simulation()->getTkp1());
            subprod(*ID, *(static_pointer_cast<LagrangianRheonomousR>(rel)->hDot()), *Yp, xcoord, false); // y += hDot
          }
          else
            RuntimeException::selfThrow("ZeroOrderHold::computeFreeOutput not yet implemented for SICONOS_OSNSP ");
        }
        // For the relation of type LagrangianScleronomousR
        if (relationSubType == ScleronomousR)
        {

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
          e = static_pointer_cast<FirstOrderLinearTIR>(rel)->e();
          F = static_pointer_cast<FirstOrderLinearTIR>(rel)->F();
        }
        else
        {
          e = static_pointer_cast<FirstOrderLinearR>(rel)->e();
          F = static_pointer_cast<FirstOrderLinearR>(rel)->F();
        }

        if (e)
          *Yp += *e;

        if (F)
        {
          SiconosVector& workZ = *inter->workZ();
          coord[3] = F->size(1);
          coord[5] = F->size(1);
          subprod(*F, workZ, *Yp, coord, false);
        }
      }

    }

    if (inter->getRelationType() == Lagrangian || inter->getRelationType() == NewtonEuler)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
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
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    DynamicalSystem& ds = **it;
    Type::Siconos dsType = Type::value(ds);

    // 1 - First Order Linear Systems
    if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(ds);
      SiconosVector& x = *d.x();
      // 1 - First Order Linear Time Invariant Systems
      // \Phi is already computed
      x = *d.workFree(); // x = xfree = Phi*xold
      if (level != LEVELMAX)
      {
        if (_PsiMap[(*it)->number()])
        {
          if (dsType == Type::FirstOrderLinearDS)
          {
            // we have to find the control interaction
            DynamicalSystemsGraph& DSG0 = *simulationLink->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);
            DynamicalSystemsGraph::OEIterator oei, oeiend;
            for (boost::tie(oei, oeiend) = DSG0.out_edges(DSG0.descriptor(*it)); oei != oeiend; ++oei)
            {
              if (DSG0.properties(*oei).forControl)
              {
                Relation& rel = *DSG0.bundle(*oei)->relation();
                computePsi(d, rel);
                break;
              }
            }
          }
          SiconosMatrix& Psi = *_PsiMap[(*it)->number()];
          prod(Psi, *d.r(), x, false); // x += Psi*u
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
  assert(!isnan(y));
  if (y <= 0)
    DEBUG_PRINT("ZeroOrderHold::addInteractionInIndexSet ACTIVATE.\n");
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
  assert(!isnan(y));
  if (y > 0)
    DEBUG_PRINT("ZeroOrderHold::removeInteractionInIndexSet DEACTIVATE.\n");
  return (y > 0);
}

void ZeroOrderHold::display()
{
  OneStepIntegrator::display();

  cout << "====== ZOH OSI display ======" << endl;
  DSIterator it;
  int itN;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    cout << "--------------------------------" << endl;
    itN = (*it)->number();
    cout << "--> Phi of dynamical system number " << itN << ": " << endl;
    if (_PhiMap[itN]) _PhiMap[itN]->display();
    else cout << "-> NULL" << endl;
    cout << "--> Psi of dynamical system number " << itN << ": " << endl;
    if (_PsiMap[itN]) _PsiMap[itN]->display();
    else cout << "-> NULL" << endl;
  }
  cout << "================================" << endl;
}
ZeroOrderHold* ZeroOrderHold::convert(OneStepIntegrator* osi)
{
  ZeroOrderHold* zeroOrderHold = dynamic_cast<ZeroOrderHold*>(osi);
  return zeroOrderHold;
}
