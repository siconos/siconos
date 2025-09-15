/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "Hem5OSI.hpp"

#include "BlockVector.hpp"
#include "EventDriven.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonEulerR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"
#include "SiconosAlgebraProd.hpp"  // for prod and subprod
#include "SiconosFortran.h"        // for Fortran to C api, fprobpointer ...
#include "Topology.hpp"

using namespace RELATION;

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// ===== Out of class objects and functions =====

namespace {  // Anonymous, local scope only

// global object and wrapping functions -> required for function plug-in and call in fortran
// routine.
std::shared_ptr<Hem5OSI> hem5_global_object{nullptr};

// // ===== Hidden implementation so we don't depend on hairer.h publically =====

// This function must have the same signature as argument FPROB  in HEM5
extern "C" void Hem5OSI_fprob_wrapper(int* IFCN, int* NQ, int* NV, int* NU, int* NL, int* LDG,
                                      int* LDF, int* LDA, int* NBLK, int* NMRC, int* NPGP,
                                      int* NPFL, int* INDGR, int* INDGC, int* INDFLR,
                                      int* INDFLC, double* time, double* q, double* v,
                                      double* u, double* xl, double* G, double* GQ, double* F,
                                      double* GQQ, double* GT, double* FL, double* QDOT,
                                      double* UDOT, double* AM) {
  return hem5_global_object->_impl->fprob(IFCN, NQ, NV, NU, NL, LDG, LDF, LDA, NBLK, NMRC,
                                          NPGP, NPFL, INDGR, INDGC, INDFLR, INDFLC, time, q, v,
                                          u, xl, G, GQ, F, GQQ, GT, FL, QDOT, UDOT, AM);
}

// This function must have the same signature as argument SOLOUT in HEM5
extern "C" void Hem5OSI_solout_wrapper(int* MODE, int* NSTEP, int* NQ, int* NV, int* NU,
                                       int* NL, int* LDG, int* LDF, int* LDA, int* LRDO,
                                       int* LIDO, siconos::fortran::hairer::fprobpointer FPROB,
                                       double* q, double* v, double* u, double* DOWK,
                                       int* IDOWK) {
  return hem5_global_object->_impl->solout(MODE, NSTEP, NQ, NV, NU, NL, LDG, LDF, LDA, LRDO,
                                           LIDO, FPROB, q, v, u, DOWK, IDOWK);
}

}  // namespace

// ===== Main class implementation ====

Hem5OSI::Hem5OSI()
    : OneStepIntegrator(OSI::HEM5OSI), _idid(0) {
#if !defined(HAS_FORTRAN)
  THROW_EXCEPTION(
      "siconos::integrators::hem5OSI: the fortran interface is not active. You can not "
      "create and use Hem5 Solver. Try to configure and reinstall siconos with "
      "WITH_FORTRAN=ON");

#endif
  _steps = 1;
  _sizeMem = 2;
  // Set levels. This may depend on the nonsmooth law and will be updated during
  // initializeWorkVectorsForInteraction(...) call.
  _levelMinForOutput = 0;
  _levelMaxForOutput = 2;
  _levelMinForInput = 1;
  _levelMaxForInput = 2;
}

void Hem5OSI::setTol(int newItol, std::vector<double>&& newRtol,
                     std::vector<double>&& newAtol) {
  _intData[4] = newItol;
  // ITOL  indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //           dimension NQ + NV + NU (ITOL=1)
  rtol = std::move(newRtol);
  atol = std::move(newAtol);
}
void Hem5OSI::setTol(int newItol, double newRtol, double newAtol) {
  _intData[4] = newItol;
  // ITOL  indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //           dimension NQ + NV + NU (ITOL=1)
  rtol[0] = newRtol;  // rtol
  atol[0] = newRtol;  // atol
}

void Hem5OSI::setMaxStepSize(double _maxStep) { rwork[5] = _maxStep; }

void Hem5OSI::setMaxNstep(int _maxNumberSteps) {
  iwork[11] = _maxNumberSteps;
}

void Hem5OSI::updateIntData() {
  //   int parameters for HEM5 are saved in vector intData.

  // 1 - _intData[0] NQ size of the position vector q
  _intData[0] = _qWork->size();

  // 2 - _intData[1] NV size of the position vector v
  _intData[1] = _vWork->size();

  // 3 - _intData[2] NU size of the external dynamic vector u
  _intData[2] = 0;

  // 4 -  _intData[3] NL size of the Lagrange multiplier vector lambda
  _intData[3] = numberOfConstraints();

  // 3 - Itol, itask, iopt
  _intData[4] = 0;  // ITOL indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //  dimension NQ + NV + NU (ITOL=1)
  _intData[5] = 0;  // IOUT selects the dense output formula

  // this computation has to be redone every time _indData[3] is recompyuted.

  // IWK(14)  MODE (=0: FULL LINEAR ALGEBRA WITH DEC, =1: IDEM WITH FL,
  //                      =2: FULL LINEAR ALGEBRA WITH DGETRF, =3: FL
  //                      =4: SPARSE, =5: IDEM WITH FL)
  int MODE = 0;
  _intData[8] = MODE;
  int NZA = 0;
  int LL = 0;
  int IS = 0;   // size of IMEM common work space arrays for MA28PACK
  int IXS = 0;  // size of XMEM common work space arrays for MA28PACK

  int LDG = 0;  // LDG : leading dimension of the Jacabian of constraints (G) (or non zeroa
                // elements in sparse case)
  int LDF = 0;  // LDF : leading dimension of the L or FL (L)

  int NMRC = (int)_intData[1];  // NMRC : size of a block of M
  int NBLK = 1;                 // NBLK : number of block of M

  if (MODE <= 3) {
    LL = 8 * ((int)_intData[1] * (int)_intData[3]) +
         4 * ((int)_intData[1] + (int)_intData[3]) * ((int)_intData[1] + (int)_intData[3]);
    LDG = _intData[3];
    LDF = _intData[3];
    NZA = LDG + std::max(LDG, LDF) + NMRC * NMRC * NBLK;
    IS = 0;   // Sparse solver MA28 is not called
    IXS = 0;  // Sparse solver MA28 is not called
  }
  if (MODE > 3) {
    THROW_EXCEPTION("Hem5OSI::updateIntData(), MODE >3 Sparse case not implemented ...");
  }

  // 5 - LWK length of real array rwork
  _intData[6] = 19 + 27 * (int)_intData[0] + 28 * (int)_intData[1] + 27 * (int)_intData[2] +
                5 * ((int)_intData[1] + (int)_intData[3]) + 4 * NZA + 2 * IXS + LL;

  // 6 - LIWK length of int array iwork
  _intData[7] =
      95 + 2 * ((int)_intData[1] + (int)_intData[3]) + 2 * IS + 12 * LDG + 4 * LDF + 4 * NZA;
  _intData[7] *= 2;
}

void Hem5OSI::updateData() {
  // Used to update some data (iwork ...) when _intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values are set.

  unsigned int sizeTol = _intData[0];  // size of rtol, atol ...
  // If itol (_intData[4]) = 0 => scalar else, vector of size neq (_intData[0]).
  //  if(_intData[0]==1) sizeTol = 1;
  //  else sizeTol = _intData[0];

  rtol.resize(sizeTol);  // rtol, relative tolerance

  atol.resize(sizeTol, 0.);  // atol, absolute tolerance

  iwork.resize(_intData[7], 0);

  rwork.resize(_intData[6], 0.);
}

void Hem5OSI::fillqWork(int* NQ, double* q) {
  unsigned int sizeQ = (unsigned int)(*NQ);
  for (unsigned int i = 0; i < sizeQ; ++i) (*_qWork)(i) = q[i];
}

void Hem5OSI::fillvWork(int* NV, double* v) {
  unsigned int sizeV = (unsigned int)(*NV);
  for (unsigned int i = 0; i < sizeV; ++i) (*_vWork)(i) = v[i];
}

void Hem5OSI::computeRhs(double t) {
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeRhs(t);
  }
}

void Hem5OSI::computeJacobianRhs(double t) {
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeJacobianRhsx(t);
  }
}
void Hem5OSI::Hem5OSI_impl::fprob(
    int* IFCN, int* NQ, int* NV, int* NU, int* NL, int* LDG, int* LDF, int* LDA, int* NBLK,
    int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR, int* INDFLC,
    double* time, double* q, double* v, double* u, double* xl, double* G, double* GQ,
    double* F, double* GQQ, double* GT, double* FL, double* QDOT, double* UDOT, double* AM) {
  DEBUG_PRINTF(
      "siconos::integrators::Hem5OSI::fprob(int* IFCN,...) with IFCN = "
      "%i \n",
      (int)*IFCN);
  DEBUG_PRINTF("NQ = %i\t NV = %i \t NU = %i, NL = %i \n", (int)*NQ, (int)*NV, (int)*NU,
               (int)*NL);
  DEBUG_PRINTF("LDG = %i\t LDF = %i \t LDA = %i \n", (int)*LDG, (int)*LDF, (int)*LDA);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  hem5osi->fillqWork(NQ, q);
  hem5osi->fillvWork(NV, v);

  double t = *time;

  auto dsGraph = hem5osi->_dynamicalSystemsGraph;

  int ifcn = (int)(*IFCN);

  if ((ifcn == 1) || (ifcn >= 7))  // compute Mass AM
  {
    unsigned int pos = 0;
    for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi) {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if (Type::value(*ds) == Type::LagrangianDS ||
          Type::value(*ds) == Type::LagrangianLinearTIDS) {
        LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds);
        if (lds.mass()) {
          lds.computeMass();
          for (unsigned int ii = pos; ii < ((unsigned int)(*NV) + pos); ii++) {
            for (unsigned int jj = pos; jj < ((unsigned int)(*NV) + pos); jj++) {
              AM[ii + jj * (int)(*NV)] = lds.mass()->getValue(ii, jj);
            }
          }
        } else {
          for (unsigned int ii = pos; ii < ((unsigned int)(*NV) + pos); ii++) {
            for (unsigned int jj = pos; jj < ((unsigned int)(*NV) + pos); jj++) {
              if (ii == jj)
                AM[ii + jj * (int)(*NV)] = 1.;
              else
                AM[ii + jj * (int)(*NV)] = 0.;
            }
          }
        }
        pos += lds.dimension();
      } else {
        THROW_EXCEPTION("Hem5OSI::fprob(), Only integration of Lagrangian DS is allowed");
      }
      DEBUG_EXPR(for (int kk = 0; kk < (int)(*NV) * (int)(*NV);
                      kk++) { std::cout << AM[kk] << std::endl; });
    }
  }
  if ((ifcn == 1) || (ifcn == 5) || (ifcn == 7) || (ifcn == 8))  // compute F
  {
    for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi) {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if (Type::value(*ds) == Type::LagrangianDS ||
          Type::value(*ds) == Type::LagrangianLinearTIDS) {
        LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds);
        hem5osi->fillqWork(NQ, q);
        hem5osi->fillvWork(NV, v);
        lds.computeForces((double)*time, lds.q(), lds.velocity());
      } else if (Type::value(*ds) == Type::NewtonEulerDS) {
        THROW_EXCEPTION(
            "Hem5OSI::fprob(), Integration of Newton Euler DS not yet implemented.");
      } else {
        THROW_EXCEPTION("Hem5OSI::fprob(), Only integration of Lagrangian DS is allowed");
      }
    }
    for (unsigned int ii = 0; ii < (unsigned int)(*NV); ii++) {
      F[ii] = hem5osi->_forcesWork->getValue(ii);
    }
  }
  if (ifcn == 4)  // compute G (constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2 =
        hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    assert(indexSet2);
    for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
      SP::Interaction inter = indexSet2->bundle(*ui);
      inter->computeOutput(t, 0);
      assert(0);
    }
  }

  if ((ifcn == 6) || (ifcn >= 10))  // compute GP ( Jacobian of the constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2 =
        hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
      SP::Interaction inter = indexSet2->bundle(*ui);
      inter->relation()->computeJach(t, *inter);
      assert(0);
    }
  }

  if ((ifcn == 5) || (ifcn == 7))  // compute GPP ( Hessian of the constraints)
  {
    // THROW_EXCEPTION("Hem5OSI::fprob(), G_qq is not available");
    std::cout << "Hem5OSI::fprob(), G_qq is not available " << std::endl;
  }

  if ((ifcn == 3) || (ifcn == 6) ||
      (ifcn >= 10))  // compute GT (partial time derivative of the constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2 =
        hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
      auto inter = indexSet2->bundle(*ui);
      inter->relation()->computeJach(t, *inter);
      assert(0);
    }
  }

  if (ifcn == 0)  // compute UDOT
  {
    for (int ii = 0; ii < (int)*NU; ii++) {
      assert(0);
    }
  }

  if ((ifcn == 1) || (ifcn == 2) || (ifcn == 10))  // compute QDOT
  {
    unsigned int pos = 0;
    for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi) {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if (Type::value(*ds) == Type::LagrangianDS ||
          Type::value(*ds) == Type::LagrangianLinearTIDS) {
        LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds);
        unsigned int dim = lds.dimension();
        for (unsigned int i = 0; i < dim; i++) {
          QDOT[i + pos] = v[i + pos];
        }
        pos += dim;
      } else if (Type::value(*ds) == Type::NewtonEulerDS) {
        THROW_EXCEPTION(
            "Hem5OSI::fprob(), Integration of Newton Euler DS not yet implemented.");
      } else {
        THROW_EXCEPTION("Hem5OSI::fprob(), Only integration of Mechanical DS is allowed");
      }
    }
    DEBUG_EXPR(for (int kk = 0; kk < (int)(*NV); kk++) { std::cout << QDOT[kk] << "\n"; });
  }

  DEBUG_PRINTF("END : Hem5OSI::fprob(integer* IFCN,...) with IFCN = %i \n \n", (int)*IFCN);
}

// void siconos::integrators::Hem5OSI::g(int* nEq,
// double*  time, double* x,
// int* ng, double* gOut)
// {
//   std::static_pointer_cast<EventDriven>(_simulation)->computeg(shared_from_this(), nEq,
//   time, x, ng, gOut);
// }

// void Hem5OSI::jacobianfx(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml,
// integer* mu,  doublereal* jacob, integer* nrowpd)
// {
//   std::static_pointer_cast<EventDriven>(_simulation)->computeJacobianfx(shared_from_this(),
//   sizeOfX, time, x, jacob);
// }
void Hem5OSI::initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) {
  // Get work buffers from the graph
  VectorOfVectors& ds_work_vectors = *_initializeDSWorkVectors(ds);

  Type::Siconos dsType = Type::value(*ds);

  ds->initRhs(t);  // This will create p[2] and other required vectors/buffers

  if (!_qWork) _qWork.reset(new BlockVector());
  if (!_vWork) _vWork.reset(new BlockVector());
  if (!_aWork) _aWork.reset(new BlockVector());
  if (!_uWork) _uWork.reset(new BlockVector());
  if (!_lambdaWork) _lambdaWork.reset(new BlockVector());
  if (!_forcesWork) _forcesWork.reset(new BlockVector());

  if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS) {
    LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds);
    lds.init_inverse_mass();  // invMass required to update post-impact velocity

    _qWork->insertPtr(lds.q());
    _vWork->insertPtr(lds.velocity());
    _aWork->insertPtr(lds.acceleration());
    _forcesWork->insertPtr(lds.forces());
    ds_work_vectors.resize(Hem5OSI::WORK_LENGTH);
    ds_work_vectors[Hem5OSI::FREE].reset(new SiconosVector(lds.dimension()));

  } else {
    THROW_EXCEPTION("Hem5OSI::initialize(), Only integration of Lagrangian DS is allowed");
  }

  ds->swapInMemory();
}

void Hem5OSI::initializeWorkVectorsForInteraction(Interaction& inter,
                                                  InteractionProperties& interProp,
                                                  DynamicalSystemsGraph& DSG) {
  SP::DynamicalSystem ds1 = interProp.source;
  SP::DynamicalSystem ds2 = interProp.target;

  if (!interProp.workVectors) {
    interProp.workVectors.reset(new VectorOfVectors);
    interProp.workVectors->resize(Hem5OSI::WORK_INTERACTION_LENGTH);
  }

  if (!interProp.workBlockVectors) {
    interProp.workBlockVectors.reset(new VectorOfBlockVectors);
    interProp.workBlockVectors->resize(Hem5OSI::BLOCK_WORK_LENGTH);
  }

  VectorOfVectors& inter_work = *interProp.workVectors;
  VectorOfBlockVectors& inter_work_block = *interProp.workBlockVectors;

  Relation& relation = *inter.relation();
  RELATION::TYPES relationType = relation.getType();

  inter_work[Hem5OSI::OSNSP_RHS].reset(new SiconosVector(inter.dimension()));

  NonSmoothLaw& nslaw = *inter.nonSmoothLaw();
  Type::Siconos nslType = Type::value(nslaw);

  if (nslType == Type::NewtonImpactNSL || nslType == Type::MultipleImpactNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  } else if (nslType == Type::NewtonImpactFrictionNSL) {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    THROW_EXCEPTION(
        "HEM5OSI::initializeWorkVectorsForInteraction  not yet implemented for nonsmooth law "
        "of type NewtonImpactFrictionNSL");
  } else
    THROW_EXCEPTION(
        "HEM5OSI::initializeWorkVectorsForInteraction not yet implemented  for nonsmooth of "
        "type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  inter.initializeMemory(_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2)))) {
    THROW_EXCEPTION(
        "siconos::integrators::Hem5OSI::initializeWorkVectorsForInteraction. The "
        "implementation is not correct for two different OSI for one interaction");
  }

  auto& workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  auto& DSlink = inter.linkToDSVariables();

  if (relationType == Lagrangian) {
    inter_work_block[Hem5OSI::xfree].reset(new BlockVector());
    inter_work_block[Hem5OSI::xfree]->insertPtr(workVds1[Hem5OSI::FREE]);
    LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds1);
    DSlink[LagrangianR::q2].reset(new BlockVector());
    DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
  }
  // else if (relationType == NewtonEuler)
  // {
  //   inter_work_block[::xfree].reset(new BlockVector());
  //   inter_work_block[::xfree]->insertPtr(workVds1[Hem5OSI::FREE]);
  // }

  if (ds1 != ds2) {
    VectorOfVectors& workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == Lagrangian) {
      inter_work_block[Hem5OSI::xfree]->insertPtr(workVds2[Hem5OSI::FREE]);
      LagrangianDS& lds = *std::static_pointer_cast<LagrangianDS>(ds2);
      DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
    }
    // else if (relationType == NewtonEuler)
    // {
    //   inter_work_block[::xfree]->insertPtr(workVds2[Hem5OSI::FREE]);
    // }
  }
}

void Hem5OSI::initialize() {
  DEBUG_PRINT("Hem5OSI::initialize(Model& m)\n");

  OneStepIntegrator::initialize();

  _impl =
      std::make_shared<Hem5OSI_impl>(std::static_pointer_cast<Hem5OSI>(shared_from_this()));

  // siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  // auto indexSet0
  //   = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  // assert(indexSet0);
  // for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  // {
  //   auto inter = indexSet0->bundle(*ui);
  //   _lambdaWork->insertPtr(inter->lambda(0));
  // }
}

void Hem5OSI::Hem5OSI_impl::solout(
    int* MODE, int* NSTEP, int* NQ, int* NV, int* NU, int* NL, int* LDG, int* LDF, int* LDA,
    int* LRDO, int* LIDO, siconos::fortran::hairer::fprobpointer FPROB, double* q, double* v,
    double* u, double* DOWK, int* IDOWK)

{}

unsigned int Hem5OSI::numberOfConstraints() {
  DEBUG_PRINT("Hem5OSI::updateConstraints() \n");
  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet2 =
      _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
  assert(indexSet2);
  SP::SiconosVector y;
  unsigned int n = 0;
  for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
    SP::Interaction inter = indexSet2->bundle(*ui);
    n++;
  }
  return n;
}

void Hem5OSI::integrate(double& tinit, double& tend, double& tout, int& idid) {
  DEBUG_PRINT(
      "Hem5OSI::integrate(double& tinit, double& tend, double& tout, int& idid) with \n");
  DEBUG_PRINTF("tinit = %f, tend= %f, tout = %f, idid = %i\n", tinit, tend, tout, idid);

  double tend_DR = tend;    // next point where output is desired (different from t!)
  double tinit_DR = tinit;  // current (starting) time

  // === Pointers to function ===
  //  --> definition and initialisation thanks to wrapper:
  hem5_global_object = std::static_pointer_cast<Hem5OSI>(
      shared_from_this());  // Warning: global object must be initialized to current one before
                            // pointers to function initialisation.

  // === HEM5 CALL ===

  updateIntData();
  if (!_qtmp) {
    _qtmp.reset(new SiconosVector(_qWork->size()));
  } else
    _qtmp->resize((int)_intData[0], true);

  DEBUG_PRINTF("Hem5OSI::integrate() _intData[0] (NQ) = %i \n", _intData[0]);

  if (!_vtmp) {
    _vtmp.reset(new SiconosVector(_vWork->size()));
  } else
    _vtmp->resize((int)_intData[1], true);

  _utmp.reset(new SiconosVector(1));
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[2] (NU) = %i \n", _intData[2]);

  if (!_atmp) {
    _atmp.reset(new SiconosVector(_vWork->size()));
  } else
    _atmp->resize((int)_intData[1], true);

  if (!_lambdatmp) {
    _lambdatmp.reset(new SiconosVector(_intData[3], 0.0));
  } else
    _lambdatmp->resize((int)_intData[3], true);
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[3] (NL) = %i \n", _intData[3]);

  DEBUG_PRINTF("Hem5OSI::integrate() _intData[6] (LWK) = %i \n", _intData[6]);
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[7] (LIWK) = %i \n", _intData[7]);

  Hem5OSI::updateData();

  rwork[0] =
      siconos::internal::MACHINE_PREC;  // WK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.

  rwork[1] = 0.0;  // WK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
  //         DEFAULT 0.85D0.
  rwork[2] = 0.0;  // WK(3), WK(4)   PARAMETERS FOR STEP SIZE SELECTION
  rwork[3] = 0.0;  //                THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
  //                WK(3) <= HNEW/HOLD <= WK(4).
  //                DEFAULT VALUES: WK(3)=0.2D0, WK(4)=10.D0
  rwork[5] = 0.0;  // WK(6)   MAXIMAL STEP SIZE, DEFAULT TEND-T.

  rwork[6] = 0.0;  // WK(7) = BETA, DEFAULT 0.D0
  rwork[7] = 0.0;  // WK(8) = ALPHA, DEFAULT 1/5

  iwork[10] = 0;  // IWK(11)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
  //          THE DEFAULT VALUE (FOR IWK(11)=0) IS 100000.
  iwork[11] = 0;  // IWK(12)  SWITCH FOR A PROJECTION TO ENSURE CONSISTENT INITIAL VALUE
  //          FOR IWK(12)=1 AN INITIAL PROJECTION IS PERFORMED.
  //          NO PROJECTION IS DONE IF IWK(12)=0.
  //          THE DEFAULT VALUE FOR IWK(12) IS 0.

  iwork[12] = 0;  // IWK(13)  FOR IWK(13).GT.0 IT IS THE NUMBER OF STEPS BETWEEN
  //          TWO PROJECTIONS ON THE MANIFOLD  DEFINED BY 0 = g(q,t).
  //          FOR IWK(13).LE.0 NO PROECTION IS PERFORMED.
  //          THE DEFAULT VALUE FOR IWK(13) IS 0.

  iwork[13] =
      _intData[8];  // IWK(14)  MODE (=0: FULL LINEAR ALGEBRA WITH DEC, =1: IDEM WITH FL,
  //                =2: FULL LINEAR ALGEBRA WITH DGETRF, =3: FL
  //                =4: SPARSE, =5: IDEM WITH FL)

  iwork[14] = 1;  // IWK(15)  IACC (=1: COMPUTE THE ACCELERATION)

  iwork[15] = 1;  // IWK(16)  IGIIN (=1: COMPUTE NUMERICALLY GII)

  // C    IWK(21->29)  IPAR
  // C    IPAR(1) = IWK(21) = NMRC (SIZE OF A BLOCK OF AM)
  // C    IPAR(2) = IWK(22) = NBLK (NUMBER OF BLOCK OF AM)
  // C    IPAR(3) = IWK(23) = NPGP (0 IF GP AS THE SAME PATTERN AS PREVIOUS CALL)
  // C    IPAR(4) = IWK(24) = NPFL (0 IF FL AS THE SAME PATTERN AS PREVIOUS CALL)
  // C    IPAR(5) = IWK(25) = IS (SIZE OF INTEGER WORK SPACE FOR MA28 (MIN 13*NM))
  // C    IPAR(6) = IWK(26) = IXS (SIZE OF REAL WORK SPACE FOR MA28 (MIN NM+4*NZA))
  // C    IPAR(7) = IWK(27) = PREVL
  // C    IPAR(8) = IWK(28) = IO
  // C    IPAR(9) = FLAG TO INDICATE IF UMDFAC HAS BEEN CALLED AT LEAST ONCE

  DEBUG_EXPR(iwork[26] = 2; printf("\n"));

  // Set atol and rtol values ...
  rtol[0] = HEM5_RTOL_DEFAULT;  // rtol
  atol[0] = HEM5_ATOL_DEFAULT;  // atol

  *_qtmp = *_qWork;  // Copy into a continuous memory chuck
  *_vtmp = *_vWork;  // Copy into a continuous memory chuck
  //*_utmp = *_uWork; // Copy into a continuous memory chuck
  *_atmp = *_aWork;  // Copy into a continuous memory chuck

  DEBUG_EXPR(_qtmp->display(););
  DEBUG_EXPR(_vtmp->display(););
  DEBUG_EXPR(_atmp->display(););

  //*_lambdatmp = *_lambdaWork; // Copy into a continuous memory chuck

  assert(_qtmp);
  assert(_vtmp);
  assert(_utmp);
  assert(_atmp);
  assert(_lambdatmp);
  assert(_intData[7]);

  // Management of vectors of Size 0
  double* pointerToU;
  if (_intData[2] == 0)
    pointerToU = nullptr;
  else
    pointerToU = &(*_utmp)(0);

  double* pointerToXL;
  if (_intData[3] == 0)
    pointerToXL = nullptr;
  else
    pointerToXL = &(*_lambdatmp)(0);
  // call HEM5 to integrate dynamical equation
  siconos::fortran::hairer::hem5(
      &(_intData[0]), &(_intData[1]), &(_intData[2]), &(_intData[3]), &Hem5OSI_fprob_wrapper,
      &tinit_DR, &(*_qtmp)(0), &(*_vtmp)(0), pointerToU, &(*_atmp)(0), pointerToXL, &tend_DR,
      &_timeStep, &rtol.front(), &atol.front(), &(_intData[4]), &Hem5OSI_solout_wrapper,
      &(_intData[5]), &rwork.front(), &(_intData[6]), &iwork.front(), &(_intData[7]), &_idid);

  // THROW_EXCEPTION(
  //     "Hem5, Fortran Language is not enabled in siconos kernel. Compile with fortran if you
  //     " "need Hem5");
  // === Post ===
  if (_idid < 0)  // if istate < 0 => HEM2 failed
  {
    std::cout << "Hem5OSI::integrate(...) failed - idid = " << _idid
              << "\n";
    std::cout << " -1 means input is not consistent" << "\n";
    std::cout << " -2 means larger NMAX needed." << "\n";
    std::cout << " -3 means step size becomes too small." << "\n";
    std::cout << " -4 means matrix is singular" << "\n";
    std::cout << " -5 means initial projection: no convergence" << "\n";
    THROW_EXCEPTION("Hem5OSI::integrate(), integration failed");
  }

  DEBUG_EXPR_WE(std::cout << "HEM5 Statitics : " << std::endl;
                std::cout << "NSTEP = " << iwork[30] << std::endl;
                std::cout << "NACCPT = " << iwork[31] << std::endl;
                std::cout << "NREJCT = " << iwork[32] << std::endl;
                std::cout << "NFCN = " << iwork[33] << std::endl;
                std::cout << "NDEC = " << iwork[34] << std::endl;
                std::cout << "NSOL = " << iwork[35] << std::endl;);
  *_qWork = *_qtmp;
  *_vWork = *_vtmp;
  *_aWork = *_atmp;

  DEBUG_PRINTF("tend_DR = %f\n", (double)tend_DR);
  DEBUG_EXPR(_qWork->display());
  DEBUG_EXPR(_vWork->display());
  DEBUG_EXPR(_aWork->display());
  DEBUG_PRINT("\n");
  DEBUG_PRINT("\n");

  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet2 =
      _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
  assert(indexSet2);
  SP::SiconosVector y;
  unsigned int pos = 0;
  for (std::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui) {
    SP::Interaction inter = indexSet2->bundle(*ui);
    inter->lambda(2)->setValue(0, (*_lambdatmp)(pos));
    pos++;
  }

  tout = tinit_DR;  // real ouput time
  tend = tend_DR;   // necessary for next start of HEM5
}

void Hem5OSI::updateState(const unsigned int level) {
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if (level == 1)  // ie impact case: compute velocity
  {
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      SP::LagrangianDS lds =
          std::static_pointer_cast<LagrangianDS>(_dynamicalSystemsGraph->bundle(*dsi));
      lds->computePostImpactVelocity();
    }
  } else if (level == 2) {
    double time = _simulation->nextTime();
    for (std::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi) {
      if (!checkOSI(dsi)) continue;
      {
        SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
        ds->update(time);
      }
    }
  } else
    THROW_EXCEPTION("Hem5OSI::updateState(index), index is out of range. Index = " +
                    std::to_string(level));
}

struct Hem5OSI::_NSLEffectOnFreeOutput : public SiconosVisitor {
  using SiconosVisitor::visit;

  OneStepNSProblem* _osnsp;
  SP::Interaction _inter;
  InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem* p, SP::Interaction inter,
                         InteractionProperties& interProp)
      : _osnsp(p), _inter(inter), _interProp(interProp){};

  void visit(const NewtonImpactNSL& nslaw) {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    SiconosVector& osnsp_rhs = *(*_interProp.workVectors)[Hem5OSI::OSNSP_RHS];
    subscal(e, _inter->y_k(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord,
            false);  // q = q + e * q
  }

  // visit function added by Son (9/11/2010)
  void visit(const MultipleImpactNSL& nslaw) { ; }
  // note : no NewtonImpactFrictionNSL
};

void Hem5OSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter,
                                OneStepNSProblem* osnsp) {
  SP::OneStepNSProblems allOSNS = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = inter->linkToDSVariables();
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;
  SP::Interaction mainInteraction = inter;
  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix C;
  //   SP::SiconosMatrix  D;
  //   SP::SiconosMatrix  F;
  SiconosVector& osnsp_rhs =
      *(*indexSet->properties(vertex_inter).workVectors)[Hem5OSI::OSNSP_RHS];
  SP::BlockVector Xfree;

  /* V.A. 10/10/2010
   * Following the type of OSNS  we need to retrieve the velocity or the acceleration
   * This tricks is not very nice but for the moment the OSNS do not known if
   * it is in accelaration of not
   */

  // SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  if (((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
    if (relationType == Lagrangian) {
      Xfree = DSlink[Hem5OSI::xfree];
    }
    // else if  (relationType == NewtonEuler)
    // {
    //   Xfree = inter->data(::FREE);
    // }
    assert(Xfree);
    //        std::cout << "Computeqblock Xfree (Gamma)========" << std::endl;
    //       Xfree->display();
  } else if (((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp) {
    Xfree = DSlink[LagrangianR::q1];
    //        std::cout << "Computeqblock Xfree (Velocity)========" << std::endl;
    //       Xfree->display();

  } else
    THROW_EXCEPTION(" computeqBlock for Event Event-driven is wrong ");

  if (relationType == Lagrangian) {
    C = mainInteraction->relation()->C();
    if (C) {
      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      subprod(*C, *Xfree, osnsp_rhs, coord, true);
    }

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
    if (relationSubType == RheonomousR) {
      if (((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
        THROW_EXCEPTION(
            "Hem5OSI::computeFreeOutput not yet implemented for LCP at acceleration level "
            "with LagrangianRheonomousR");
      } else if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp) {
        std::static_pointer_cast<LagrangianRheonomousR>(inter->relation())
            ->computehDot(simulation()->getTkp1(), *DSlink[LagrangianR::q0],
                          *DSlink[LagrangianR::z]);
        subprod(*ID,
                *(std::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()),
                osnsp_rhs, xcoord, false);  // y += hDot
      } else
        THROW_EXCEPTION("Hem5OSI::computeFreeOutput not implemented for SICONOS_OSNSP ");
    }
    // For the relation of type LagrangianScleronomousR
    if (relationSubType == ScleronomousR) {
      if (((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp) {
        std::static_pointer_cast<LagrangianScleronomousR>(inter->relation())
            ->computedotjacqhXqdot(simulation()->getTkp1(), *inter, DSlink);
        subprod(*ID,
                *(std::static_pointer_cast<LagrangianScleronomousR>(inter->relation())
                      ->dotjacqhXqdot()),
                osnsp_rhs, xcoord, false);  // y += NonLinearPart
      }
    }
  } else
    THROW_EXCEPTION("Hem5OSI::computeFreeOutput not yet implemented for Relation of type " +
                    std::to_string(relationType));
  if (((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp) {
    if (inter->relation()->getType() == Lagrangian ||
        inter->relation()->getType() == NewtonEuler) {
      SP::SiconosVisitor nslEffectOnFreeOutput(
          new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }
}
void Hem5OSI::display()  {
  OneStepIntegrator::display();
  std::cout << " --- > Hem5OSI specific values: " << std::endl;
  std::cout << "Number of equations: " << _intData[0] << std::endl;
  std::cout << "Number of constraints: " << _intData[1] << std::endl;
  std::cout << "itol, itask, istate, iopt, lrw, liw, jt: (for details on what are these "
               "variables see opkdmain.f)"
            << std::endl;
  std::cout << _intData[2] << ", " << _intData[3] << ", " << _intData[4] << ", " << _intData[5]
            << ", " << _intData[6] << ", " << _intData[7] << ", " << _intData[8] << std::endl;
  std::cout << "====================================" << std::endl;
}
