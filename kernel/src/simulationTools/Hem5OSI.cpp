/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <hairer.h>
#include "EventDriven.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "BlockVector.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Model.hpp"
#include "Topology.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianScleronomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonEulerR.hpp"
#include "OneStepNSProblem.hpp"

using namespace RELATION;

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

// initial step size guess (typical value 1e-3)
#define INITIAL_GUESS_TS 1.e-3

// ===== Hidden implementation so we don't depend on hairer.h publically =====

class Hem5OSI_impl
{
public:
  Hem5OSI_impl(Hem5OSI* h) : hem5osi(h) {}
  Hem5OSI *hem5osi;
  fprobfunction fprob;
  soloutfunction solout;
};

// ===== Out of class objects and functions =====

// global object and wrapping functions -> required for function plug-in and call in fortran routine.
SP::Hem5OSI hem5_global_object;

// This first function must have the same signature as argument FPROB  in HEM5
extern "C" fprobfunction Hem5OSI_fprob_wrapper;

void Hem5OSI_fprob_wrapper(integer* IFCN,
                           integer* NQ,
                           integer* NV,
                           integer* NU,
                           integer* NL,
                           integer* LDG, integer* LDF, integer* LDA,
                           integer* NBLK, integer* NMRC,
                           integer* NPGP, integer* NPFL,
                           integer* INDGR, integer* INDGC, integer * INDFLR, integer * INDFLC,
                           doublereal* time,
                           doublereal* q, doublereal* v, doublereal* u,  doublereal* xl,
                           doublereal* G, doublereal* GQ, doublereal * F,
                           doublereal* GQQ, doublereal* GT, doublereal * FL,
                           doublereal* QDOT, doublereal* UDOT, doublereal * AM)
{
  return hem5_global_object->_impl->fprob(IFCN,
                                   NQ,
                                   NV,
                                   NU,
                                   NL,
                                   LDG,  LDF,  LDA,
                                   NBLK,  NMRC,
                                   NPGP,  NPFL,
                                   INDGR,  INDGC, INDFLR, INDFLC,
                                   time,
                                   q,  v,  u,   xl,
                                   G,  GQ, F,
                                   GQQ,  GT, FL,
                                   QDOT,  UDOT, AM);
}

// This first function must have the same signature as argument SOLOUT in HEM5
extern "C" soloutfunction Hem5OSI_solout_wrapper;
void Hem5OSI_solout_wrapper(integer* MODE,
                            integer* NSTEP,
                            integer* NQ,
                            integer* NV,
                            integer* NU,
                            integer* NL,
                            integer* LDG, integer* LDF, integer* LDA,
                            integer* LRDO, integer* LIDO,
                            fprobpointer FPROB,
                            doublereal* q, doublereal* v, doublereal* u,
                            doublereal *DOWK, integer* IDOWK)
{
  return hem5_global_object->_impl->solout(MODE,
                                    NSTEP,
                                    NQ,
                                    NV,
                                    NU,
                                    NL,
                                    LDG,  LDF, LDA,
                                    LRDO, LIDO,
                                    FPROB,
                                    q, v,  u,
                                    DOWK, IDOWK);
}

// ===== Main class implementation ====

Hem5OSI::Hem5OSI():
  OneStepIntegrator(OSI::HEM5OSI), _idid(0)
  , _impl(std11::make_shared<Hem5OSI_impl>(this))
{
  _steps=1;
  _intData.resize(9);
  for(int i = 0; i < 9; i++) _intData[i] = 0;
  _sizeMem = 2;
  _timeStep = INITIAL_GUESS_TS;
  // Set levels. This may depend on the nonsmooth law and will be updated during fillDSLinks(...) call.
  _levelMinForOutput=0;
  _levelMaxForOutput=2;
  _levelMinForInput=1;
  _levelMaxForInput=2;
}

void Hem5OSI::setTol(integer newItol, SA::doublereal newRtol, SA::doublereal newAtol)
{
  _intData[4] = newItol; // ITOL  indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //           dimension NQ + NV + NU (ITOL=1)
  rtol = newRtol;
  atol = newAtol;
}
void Hem5OSI::setTol(integer newItol, doublereal newRtol, doublereal newAtol)
{
  _intData[4] = newItol; // ITOL  indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //           dimension NQ + NV + NU (ITOL=1)
  rtol[0] = newRtol; // rtol
  atol[0] = newRtol;  // atol
}

void Hem5OSI::setMaxStepSize(doublereal _maxStep)
{
  rwork[5] = _maxStep;
}

void Hem5OSI::setMaxNstep(integer _maxNumberSteps)
{
  iwork[11] = _maxNumberSteps;
}

void Hem5OSI::updateIntData()
{
  //   Integer parameters for HEM5 are saved in vector intData.

  // 1 - _intData[0] NQ size of the position vector q
  _intData[0] = _qWork->size();

  // 2 - _intData[1] NV size of the position vector v
  _intData[1] = _vWork->size();

  // 3 - _intData[2] NU size of the external dynamic vector u
  _intData[2] = 0;

  // 4 -  _intData[3] NL size of the Lagrange multiplier vector lambda
  _intData[3] = numberOfConstraints();

  // 3 - Itol, itask, iopt
  _intData[4] = 0; // ITOL indicates whether RTOL and ATOL are scalar (ITOL=0), or array of
  //  dimension NQ + NV + NU (ITOL=1)
  _intData[5] = 0; // IOUT selects the dense output formula

  // this computation has to be redone every time _indData[3] is recompyuted.


  // IWK(14)  MODE (=0: FULL LINEAR ALGEBRA WITH DEC, =1: IDEM WITH FL,
  //                      =2: FULL LINEAR ALGEBRA WITH DGETRF, =3: FL
  //                      =4: SPARSE, =5: IDEM WITH FL)
  int MODE = 0;
  _intData[8] = MODE;
  int NZA = 0;
  int LL =0;
  int IS = 0;  // size of IMEM common work space arrays for MA28PACK
  int IXS = 0; // size of XMEM common work space arrays for MA28PACK

  int LDG = 0; // LDG : leading dimension of the Jacabian of constraints (G) (or non zeroa elements in sparse case)
  int LDF = 0; // LDF : leading dimension of the L or FL (L)

  int NMRC = (int)_intData[1]; // NMRC : size of a block of M
  int NBLK = 1;                // NBLK : number of block of M

  if(MODE <=3)
  {
    LL = 8 * ((int)_intData[1] * (int)_intData[3])
      + 4 * ((int)_intData[1] + (int)_intData[3])*((int)_intData[1] + (int)_intData[3]);
    LDG = _intData[3];
    LDF = _intData[3];
    NZA = LDG + std::max(LDG,LDF) + NMRC*NMRC*NBLK;
    IS  = 0; // Sparse solver MA28 is not called
    IXS = 0; // Sparse solver MA28 is not called
  }
  if(MODE >3)
  {
    RuntimeException::selfThrow("Hem5OSI::updateIntData(), MODE >3 Sparse case not implemented ...");
  }

  // 5 - LWK length of real array rwork
  _intData[6] = 19 + 27*(int)_intData[0] + 28 * (int)_intData[1] + 27 * (int)_intData[2]
    + 5*((int)_intData[1] + (int)_intData[3]) + 4*NZA + 2*IXS + LL;

  // 6 - LIWK length of integer array iwork
  _intData[7] = 95 + 2*((int)_intData[1]+(int)_intData[3]) + 2*IS + 12*LDG + 4 * LDF + 4 *NZA;
  _intData[7] *= 2;
}

void Hem5OSI::updateData()
{
  // Used to update some data (iwork ...) when _intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values are set.

  unsigned int sizeTol = _intData[0]; // size of rtol, atol ...
  // If itol (_intData[4]) = 0 => scalar else, vector of size neq (_intData[0]).
  //  if(_intData[0]==1) sizeTol = 1;
  //  else sizeTol = _intData[0];

  rtol.reset(new doublereal[sizeTol]) ;    // rtol, relative tolerance

  atol.reset(new doublereal[sizeTol]) ;    // atol, absolute tolerance
  for(unsigned int i = 0; i < sizeTol; i++)
  {
    atol[i] = 0.0;
  }

  iwork.reset(new integer[_intData[7]]);
  for(int i = 0; i < _intData[7]; i++) iwork[i] = 0;

  rwork.reset(new doublereal[_intData[6]]);
  for(int i = 0; i < _intData[6]; i++) rwork[i] = 0.0;

}

void Hem5OSI::fillqWork(integer* NQ, doublereal* q)
{
  unsigned int sizeQ = (unsigned int)(*NQ);
  for(unsigned int i = 0; i < sizeQ ; ++i)
    (*_qWork)(i) = q[i];
}

void Hem5OSI::fillvWork(integer* NV, doublereal* v)
{
  unsigned int sizeV = (unsigned int)(*NV);
  for(unsigned int i = 0; i < sizeV ; ++i)
    (*_vWork)(i) = v[i];
}

void Hem5OSI::computeRhs(double t)
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeRhs(t);
  }
}

void Hem5OSI::computeJacobianRhs(double t)
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if(!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    ds->computeJacobianRhsx(t);
  }
}

void Hem5OSI_impl::fprob(integer* IFCN,
                    integer* NQ,
                    integer* NV,
                    integer* NU,
                    integer* NL,
                    integer* LDG, integer* LDF, integer* LDA,
                    integer* NBLK, integer* NMRC,
                    integer* NPGP, integer* NPFL,
                    integer* INDGR, integer* INDGC, integer * INDFLR, integer * INDFLC,
                    doublereal* time,
                    doublereal* q, doublereal* v, doublereal* u,  doublereal* xl,
                    doublereal* G, doublereal* GQ, doublereal * F,
                    doublereal* GQQ, doublereal* GT, doublereal * FL,
                    doublereal* QDOT, doublereal* UDOT, doublereal * AM)
{
  DEBUG_PRINTF("Hem5OSI::fprob(integer* IFCN,...) with IFCN = %i \n", (int)*IFCN);
  DEBUG_PRINTF("NQ = %i\t NV = %i \t NU = %i, NL = %i \n", (int)*NQ, (int)*NV, (int)*NU, (int)*NL);
  DEBUG_PRINTF("LDG = %i\t LDF = %i \t LDA = %i \n", (int)*LDG, (int)*LDF, (int)*LDA);

  // fill in xWork vector (ie all the x of the ds of this osi) with x
  hem5osi->fillqWork(NQ, q);
  hem5osi->fillvWork(NV, v);

  double t = *time;

  SP::DynamicalSystemsGraph dsGraph =  hem5osi->_dynamicalSystemsGraph;



  int ifcn = (int)(*IFCN);

  if((ifcn == 1) || (ifcn >= 7))  // compute Mass AM
  {
    unsigned int pos=0;
    for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if(Type::value(*ds) == Type::LagrangianDS ||
         Type::value(*ds) == Type::LagrangianLinearTIDS)
	    {
	      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(ds);
	      if(lds.mass())
        {
          lds.computeMass();
          for(unsigned int ii =pos ; ii < ((unsigned int)(*NV)+pos); ii ++)
          {
            for(unsigned int jj =pos ; jj < ((unsigned int)(*NV)+pos); jj ++)
            {
              AM[ii + jj*(int)(*NV)] = lds.mass()->getValue(ii,jj) ;
            }
          }
        }
	      else
        {
          for(unsigned int ii =pos ; ii < ((unsigned int)(*NV)+pos); ii ++)
          {
            for(unsigned int jj =pos ; jj < ((unsigned int)(*NV)+pos); jj ++)
            {
              if(ii == jj)
                AM[ii + jj*(int)(*NV)] = 1.;
              else
                AM[ii + jj*(int)(*NV)] = 0.;
            }
          }
        }
	      pos += lds.dimension();
	    }
      else
	    {
	      RuntimeException::selfThrow("Hem5OSI::fprob(), Only integration of Lagrangian DS is allowed");
	    }
      DEBUG_EXPR(
        for(int kk =0 ; kk < (int)(*NV)* (int)(*NV); kk ++)
        {
          std::cout << AM[kk] << std::endl;
        }
        );
    }
  }
  if((ifcn ==1) || (ifcn == 5) || (ifcn == 7) || (ifcn==8))  // compute F
  {
    for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if(Type::value(*ds) == Type::LagrangianDS ||
         Type::value(*ds) == Type::LagrangianLinearTIDS)
	    {
	      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(ds);
	      hem5osi->fillqWork(NQ,q);
	      hem5osi->fillvWork(NV,v);
	      lds.computeForces((double)*time, lds.q(), lds.velocity());
	    }
      else if(Type::value(*ds) == Type::NewtonEulerDS)
	    {
	      RuntimeException::selfThrow("Hem5OSI::fprob(), Integration of Newton Euler DS not yet implemented.");
	    }
      else
	    {
	      RuntimeException::selfThrow("Hem5OSI::fprob(), Only integration of Lagrangian DS is allowed");
	    }
    }
    for(unsigned int ii =0 ; ii < (unsigned int)(*NV); ii ++)
    {
      F[ii] = hem5osi->_forcesWork->getValue(ii) ;
    }
  }
  if(ifcn == 4)  // compute G (constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2
      = hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    assert(indexSet2);
    for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet2->bundle(*ui);
      inter->computeOutput(t, indexSet2->properties(*ui), 0);
      assert(0);
    }

  }

  if((ifcn == 6) || (ifcn >= 10))   // compute GP ( Jacobian of the constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2 =
      hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet2->bundle(*ui);
      inter->relation()->computeJach(t, *inter, indexSet2->properties(*ui));
      assert(0);
    }
  }

  if((ifcn == 5) || (ifcn == 7))   // compute GPP ( Hessian of the constraints)
  {
    //RuntimeException::selfThrow("Hem5OSI::fprob(), G_qq is not available");
    std::cout << "Hem5OSI::fprob(), G_qq is not available " << std::endl;
  }

  if((ifcn == 3) || (ifcn == 6) || (ifcn >= 10))   // compute GT (partial time derivative of the constraints)
  {
    InteractionsGraph::VIterator ui, uiend;
    SP::InteractionsGraph indexSet2 =
      hem5osi->_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet2->bundle(*ui);
      inter->relation()->computeJach(t, *inter, indexSet2->properties(*ui));
      assert(0);
    }
  }

  if(ifcn == 0)  // compute UDOT
  {
    for(int ii = 0; ii < (int)*NU ; ii++)
    {
      assert(0);
    }
  }

  if((ifcn == 1) || (ifcn == 2) || (ifcn == 10))   // compute QDOT
  {
    unsigned int pos=0;
    for(DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
    {
      SP::DynamicalSystem ds = dsGraph->bundle(*vi);
      if(Type::value(*ds) == Type::LagrangianDS ||
         Type::value(*ds) == Type::LagrangianLinearTIDS)
	    {
	      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(ds);
	      unsigned int dim = lds.dimension();
	      for(unsigned int i =0 ; i < dim ; i++)
        {
          QDOT[i+pos] = v[i+pos];
        }
	      pos +=dim ;
	    }
      else if(Type::value(*ds) == Type::NewtonEulerDS)
	    {
	      RuntimeException::selfThrow("Hem5OSI::fprob(), Integration of Newton Euler DS not yet implemented.");
	    }
      else
	    {
	      RuntimeException::selfThrow("Hem5OSI::fprob(), Only integration of Mechanical DS is allowed");
	    }

    }
    DEBUG_EXPR(
      for(int kk =0 ; kk < (int)(*NV); kk ++)
      {
        std::cout << QDOT[kk] << std::endl;
      }
      );
  }

  DEBUG_PRINTF("END : Hem5OSI::fprob(integer* IFCN,...) with IFCN = %i \n \n", (int)*IFCN);
}
// void Hem5OSI::g(integer* nEq, doublereal*  time, doublereal* x, integer* ng, doublereal* gOut)
// {
//   std11::static_pointer_cast<EventDriven>(_simulation)->computeg(shared_from_this(), nEq, time, x, ng, gOut);
// }

// void Hem5OSI::jacobianfx(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd)
// {
//   std11::static_pointer_cast<EventDriven>(_simulation)->computeJacobianfx(shared_from_this(), sizeOfX, time, x, jacob);
// }
void Hem5OSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  VectorOfVectors& workVectors = *_initializeDSWorkVectors(ds);

  Type::Siconos dsType = Type::value(*ds);

  if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(ds);
    lds.init_inverse_mass(); // invMass required to update post-impact velocity

    _qWork->insertPtr(lds.q());
    _vWork->insertPtr(lds.velocity());
    _aWork->insertPtr(lds.acceleration());
    _forcesWork->insertPtr(lds.forces());
    workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
    workVectors[OneStepIntegrator::free].reset(new SiconosVector(lds.dimension()));

  }
  else
  {
    RuntimeException::selfThrow("Hem5OSI::initialize(), Only integration of Lagrangian DS is allowed");
  }

  ds->swapInMemory();
}


void Hem5OSI::fillDSLinks(Interaction &inter,
                            InteractionProperties& interProp,
                            DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;

  assert(interProp.DSlink);

  VectorOfVectors& workV = *interProp.workVectors;
  workV.resize(Hem5OSI::WORK_INTERACTION_LENGTH);
  workV[Hem5OSI::OSNSP_RHS].reset(new SiconosVector(inter.getSizeOfY()));
  
  VectorOfBlockVectors& DSlink = *interProp.DSlink;
  // VectorOfVectors& workVInter = *interProp.workVectors;
  // VectorOfSMatrices& workMInter = *interProp.workMatrices;

  Relation &relation =  *inter.relation();
  NonSmoothLaw & nslaw = *inter.nonSmoothLaw();
  RELATION::TYPES relationType = relation.getType();
  Type::Siconos nslType = Type::value(nslaw);

  if (nslType == Type::NewtonImpactNSL || nslType == Type::MultipleImpactNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 2 ;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
  }
  else if (nslType ==  Type::NewtonImpactFrictionNSL)
  {
    _levelMinForOutput = 0;
    _levelMaxForOutput = 4;
    _levelMinForInput = 1;
    _levelMaxForInput = 2;
    RuntimeException::selfThrow("HEM5OSI::fillDSLinks  not yet implemented for nonsmooth law of type NewtonImpactFrictionNSL");
  }
  else
    RuntimeException::selfThrow("HEM5OSI::fillDSLinks not yet implemented  for nonsmooth of type");

  // Check if interations levels (i.e. y and lambda sizes) are compliant with the current osi.
  _check_and_update_interaction_levels(inter);
  // Initialize/allocate memory buffers in interaction.
  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);

  /* allocate and set work vectors for the osi */
  if (!(checkOSI(DSG.descriptor(ds1)) && checkOSI(DSG.descriptor(ds2))))
  {
    RuntimeException::selfThrow("LsodarOSI::fillDSLinks. The implementation is not correct for two different OSI for one interaction");
  }

  VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;
  if (relationType == Lagrangian)
  {
    DSlink[LagrangianR::xfree].reset(new BlockVector());
    DSlink[LagrangianR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  }
  // else if (relationType == NewtonEuler)
  // {
  //   DSlink[NewtonEulerR::xfree].reset(new BlockVector());
  //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
  // }

  if (ds1 != ds2)
  {
    VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
    if (relationType == Lagrangian)
    {
      DSlink[LagrangianR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    }
    // else if (relationType == NewtonEuler)
    // {
    //   DSlink[NewtonEulerR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
    // }
  }
}


void Hem5OSI::initialize(Model& m)
{

  DEBUG_PRINT("Hem5OSI::initialize(Model& m)\n");

  _qWork.reset(new BlockVector());
  _vWork.reset(new BlockVector());
  _aWork.reset(new BlockVector());
  _uWork.reset(new BlockVector());
  _lambdaWork.reset(new BlockVector());
  _forcesWork.reset(new BlockVector());
  OneStepIntegrator::initialize(m);

  // InteractionsGraph::VIterator ui, uiend;
  // SP::InteractionsGraph indexSet0
  //   = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  // assert(indexSet0);
  // for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  // {
  //   SP::Interaction inter = indexSet0->bundle(*ui);
  //   _lambdaWork->insertPtr(inter->lambda(0));
  // }


}

void Hem5OSI_impl::solout(integer* MODE,
                     integer* NSTEP,
                     integer* NQ,
                     integer* NV,
                     integer* NU,
                     integer* NL,
                     integer* LDG, integer* LDF, integer* LDA,
                     integer* LRDO, integer* LIDO,
                     fprobpointer FPROB,
                     doublereal* q, doublereal* v, doublereal* u,
                     doublereal *DOWK, integer* IDOWK)

{
}

unsigned int Hem5OSI::numberOfConstraints()
{
  DEBUG_PRINT("Hem5OSI::updateConstraints() \n");
  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet2
    = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
  assert(indexSet2);
  SP::SiconosVector y;
  unsigned int n = 0;
  for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet2->bundle(*ui);
    n++;
  }
  return n;
}

void Hem5OSI::integrate(double& tinit, double& tend, double& tout, int& idid)
{

  DEBUG_PRINT("Hem5OSI::integrate(double& tinit, double& tend, double& tout, int& idid) with \n");
  DEBUG_PRINTF("tinit = %f, tend= %f, tout = %f, idid = %i\n", tinit, tend,  tout, idid);

  doublereal tend_DR = tend  ;       // next point where output is desired (different from t!)
  doublereal tinit_DR = tinit;       // current (starting) time

  // === Pointers to function ===
  //  --> definition and initialisation thanks to wrapper:
  hem5_global_object = std11::static_pointer_cast<Hem5OSI>(shared_from_this()); // Warning: global object must be initialized to current one before pointers to function initialisation.

  // function to compute the system to simulation
  fprobpointer pointerToFPROB = Hem5OSI_fprob_wrapper;

  // function to compute the system to simulation
  soloutpointer pointerToSOLOUT = Hem5OSI_solout_wrapper;

  // === HEM5 CALL ===


  updateIntData();
  if(!_qtmp)
  {
    _qtmp.reset(new SiconosVector(_qWork->size()));
  }
  else
    _qtmp->resize((int)_intData[0],true);

  DEBUG_PRINTF("Hem5OSI::integrate() _intData[0] (NQ) = %i \n",_intData[0]);

  if(!_vtmp)
  {
    _vtmp.reset(new SiconosVector(_vWork->size()));
  }
  else
    _vtmp->resize((int)_intData[1],true);


  _utmp.reset(new SiconosVector(1));
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[2] (NU) = %i \n",_intData[2]);

  if(!_atmp)
  {
    _atmp.reset(new SiconosVector(_vWork->size()));
  }
  else
    _atmp->resize((int)_intData[1],true);

  if(!_lambdatmp)
  {
    _lambdatmp.reset(new SiconosVector(_intData[3],0.0));
  }
  else
    _lambdatmp->resize((int)_intData[3],true);
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[3] (NL) = %i \n",_intData[3]);

  DEBUG_PRINTF("Hem5OSI::integrate() _intData[6] (LWK) = %i \n",_intData[6]);
  DEBUG_PRINTF("Hem5OSI::integrate() _intData[7] (LIWK) = %i \n",_intData[7]);

  Hem5OSI::updateData();

  rwork[0] = MACHINE_PREC ; // WK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.

  rwork[1] = 0.0 ;          // WK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION,
  //         DEFAULT 0.85D0.
  rwork[2] = 0.0 ; // WK(3), WK(4)   PARAMETERS FOR STEP SIZE SELECTION
  rwork[3] = 0.0 ; //                THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION
  //                WK(3) <= HNEW/HOLD <= WK(4).
  //                DEFAULT VALUES: WK(3)=0.2D0, WK(4)=10.D0
  rwork[5] = 0.0 ; // WK(6)   MAXIMAL STEP SIZE, DEFAULT TEND-T.

  rwork[6] = 0.0 ; // WK(7) = BETA, DEFAULT 0.D0
  rwork[7] = 0.0 ; // WK(8) = ALPHA, DEFAULT 1/5

  iwork[10] = 0 ; // IWK(11)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
  //          THE DEFAULT VALUE (FOR IWK(11)=0) IS 100000.
  iwork[11] = 0 ; // IWK(12)  SWITCH FOR A PROJECTION TO ENSURE CONSISTENT INITIAL VALUE
  //          FOR IWK(12)=1 AN INITIAL PROJECTION IS PERFORMED.
  //          NO PROJECTION IS DONE IF IWK(12)=0.
  //          THE DEFAULT VALUE FOR IWK(12) IS 0.

  iwork[12] = 0 ; // IWK(13)  FOR IWK(13).GT.0 IT IS THE NUMBER OF STEPS BETWEEN
  //          TWO PROJECTIONS ON THE MANIFOLD  DEFINED BY 0 = g(q,t).
  //          FOR IWK(13).LE.0 NO PROECTION IS PERFORMED.
  //          THE DEFAULT VALUE FOR IWK(13) IS 0.


  iwork[13] = _intData[8] ; // IWK(14)  MODE (=0: FULL LINEAR ALGEBRA WITH DEC, =1: IDEM WITH FL,
  //                =2: FULL LINEAR ALGEBRA WITH DGETRF, =3: FL
  //                =4: SPARSE, =5: IDEM WITH FL)

  iwork[14] = 1    ; // IWK(15)  IACC (=1: COMPUTE THE ACCELERATION)

  iwork[15] = 1    ; // IWK(16)  IGIIN (=1: COMPUTE NUMERICALLY GII)

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

  DEBUG_EXPR(iwork[26] =2; printf("\n"));

  // Set atol and rtol values ...
  rtol[0] = HEM5_RTOL_DEFAULT ; // rtol
  atol[0] = HEM5_ATOL_DEFAULT ;  // atol

  *_qtmp = *_qWork; // Copy into a continuous memory chuck
  *_vtmp = *_vWork; // Copy into a continuous memory chuck
  //*_utmp = *_uWork; // Copy into a continuous memory chuck
  *_atmp = *_aWork; // Copy into a continuous memory chuck

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
  doublereal * pointerToU;
  if(_intData[2] ==0)
    pointerToU = NULL;
  else
    pointerToU = &(*_utmp)(0);

  doublereal * pointerToXL;
  if(_intData[3] ==0)
    pointerToXL = NULL;
  else
    pointerToXL = &(*_lambdatmp)(0);

  // call HEM5 to integrate dynamical equation
  CNAME(hem5)(&(_intData[0]),
              &(_intData[1]),
              &(_intData[2]),
              &(_intData[3]),
              pointerToFPROB,
              &tinit_DR,
              &(*_qtmp)(0),
              &(*_vtmp)(0),
              pointerToU,
              &(*_atmp)(0),
              pointerToXL,
              &tend_DR,
              &_timeStep,
              rtol.get(),
              atol.get(),
              &(_intData[4]),
              pointerToSOLOUT,
              &(_intData[5]),
              rwork.get(),
              &(_intData[6]),
              iwork.get(),
              &(_intData[7]),
              &_idid);

  // === Post ===
  if(_idid < 0)  // if istate < 0 => HEM2 failed
  {
    std::cout << "Hem5OSI::integrate(...) failed - idid = " << _idid <<std::endl;
    std::cout << " -1 means input is not consistent" <<std::endl;
    std::cout << " -2 means larger NMAX needed." <<std::endl;
    std::cout << " -3 means step size becomes too small." <<std::endl;
    std::cout << " -4 means matrix is singular" <<std::endl;
    std::cout << " -5 means initial projection: no convergence" <<std::endl;
    RuntimeException::selfThrow("Hem5OSI::integrate(), integration failed");
  }

  DEBUG_EXPR_WE(std::cout << "HEM5 Statitics : " <<std::endl;
                std::cout << "NSTEP = " << iwork[30] <<std::endl;
                std::cout << "NACCPT = " << iwork[31] <<std::endl;
                std::cout << "NREJCT = " << iwork[32] <<std::endl;
                std::cout << "NFCN = " << iwork[33] <<std::endl;
                std::cout << "NDEC = " << iwork[34] <<std::endl;
                std::cout << "NSOL = " << iwork[35] <<std::endl;);
  *_qWork = *_qtmp;
  *_vWork = *_vtmp;
  *_aWork = *_atmp;

  DEBUG_PRINTF("tend_DR = %f\n", (double) tend_DR);
  DEBUG_EXPR(_qWork->display());
  DEBUG_EXPR(_vWork->display());
  DEBUG_EXPR(_aWork->display());
  DEBUG_PRINT("\n");
  DEBUG_PRINT("\n");



  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet2
    = _simulation->nonSmoothDynamicalSystem()->topology()->indexSet(2);
  assert(indexSet2);
  SP::SiconosVector y;
  unsigned int pos=0;
  for(std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet2->bundle(*ui);
    inter->lambda(2)->setValue(0,(*_lambdatmp)(pos));
    pos++;
  }

  tout  = tinit_DR; // real ouput time
  tend  = tend_DR;  // necessary for next start of HEM5

}


void Hem5OSI::updateState(const unsigned int level)
{
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  DynamicalSystemsGraph::VIterator dsi, dsend;
  if(level == 1)  // ie impact case: compute velocity
  {
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::LagrangianDS lds = std11::static_pointer_cast<LagrangianDS>(_dynamicalSystemsGraph->bundle(*dsi));
      lds->computePostImpactVelocity();
    }
  }
  else if(level == 2)
  {
    double time = _simulation->nextTime();
    for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      {
        SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
        ds->update(time);
      }
    }
  }
  else RuntimeException::selfThrow("Hem5OSI::updateState(index), index is out of range. Index = " + level);
}

struct Hem5OSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;
  InteractionProperties& _interProp;
  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter, InteractionProperties& interProp) :
    _osnsp(p), _inter(inter), _interProp(interProp)  {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    SiconosVector & osnsp_rhs = *(*_interProp.workVectors)[Hem5OSI::OSNSP_RHS];
    subscal(e, *_inter->yOld(_osnsp->inputOutputLevel()), osnsp_rhs, subCoord, false); // q = q + e * q
  }

  // visit function added by Son (9/11/2010)
  void visit(const MultipleImpactNSL& nslaw)
  {
    ;
  }
  // note : no NewtonImpactFrictionNSL
};

void Hem5OSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem * osnsp)
{
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
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
  SP::SiconosMatrix  C;
  //   SP::SiconosMatrix  D;
  //   SP::SiconosMatrix  F;
  SiconosVector& osnsp_rhs = *(*indexSet->properties(vertex_inter).workVectors)[Hem5OSI::OSNSP_RHS];
  SP::BlockVector Xfree;


  /* V.A. 10/10/2010
   * Following the type of OSNS  we need to retrieve the velocity or the acceleration
   * This tricks is not very nice but for the moment the OSNS do not known if
   * it is in accelaration of not
   */

  //SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
  {
    if(relationType == Lagrangian)
    {
      Xfree = DSlink[LagrangianR::xfree];
    }
    // else if  (relationType == NewtonEuler)
    // {
    //   Xfree = inter->data(NewtonEulerR::free);
    // }
    assert(Xfree);
    //        std::cout << "Computeqblock Xfree (Gamma)========" << std::endl;
    //       Xfree->display();
  }
  else  if(((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp)
  {

    Xfree = DSlink[LagrangianR::q1];
    //        std::cout << "Computeqblock Xfree (Velocity)========" << std::endl;
    //       Xfree->display();

  }
  else
    RuntimeException::selfThrow(" computeqBlock for Event Event-driven is wrong ");

  if(relationType == Lagrangian)
  {
    C = mainInteraction->relation()->C();
    if(C)
    {
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
    if(relationSubType == RheonomousR)
    {
      if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
	    {
	      RuntimeException::selfThrow("Hem5OSI::computeFreeOutput not yet implemented for LCP at acceleration level with LagrangianRheonomousR");
	    }
      else if(((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
	    {
	      SiconosVector q = *DSlink[LagrangianR::q0];
	      SiconosVector z = *DSlink[LagrangianR::z];
	      std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->computehDot(simulation()->getTkp1(), q, z);
	      *DSlink[LagrangianR::z] = z;
	      subprod(*ID, *(std11::static_pointer_cast<LagrangianRheonomousR>(inter->relation())->hDot()), osnsp_rhs, xcoord, false); // y += hDot
	    }
      else
        RuntimeException::selfThrow("Hem5OSI::computeFreeOutput not implemented for SICONOS_OSNSP ");
    }
    // For the relation of type LagrangianScleronomousR
    if(relationSubType == ScleronomousR)
    {
      if(((*allOSNS)[SICONOS_OSNSP_ED_SMOOTH_ACC]).get() == osnsp)
	    {
	      std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->computedotjacqhXqdot(simulation()->getTkp1(), *inter, DSlink);
	      subprod(*ID, *(std11::static_pointer_cast<LagrangianScleronomousR>(inter->relation())->dotjacqhXqdot()), osnsp_rhs, xcoord, false); // y += NonLinearPart
	    }
    }
  }
  else
    RuntimeException::selfThrow("Hem5OSI::computeFreeOutput not yet implemented for Relation of type " + relationType);
  if(((*allOSNS)[SICONOS_OSNSP_ED_IMPACT]).get() == osnsp)
  {
    if(inter->relation()->getType() == Lagrangian || inter->relation()->getType() == NewtonEuler)
    {
      SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter, indexSet->properties(vertex_inter)));
      inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
    }
  }

}
void Hem5OSI::display()
{
  OneStepIntegrator::display();
  std::cout << " --- > Hem5OSI specific values: " <<std::endl;
  std::cout << "Number of equations: " << _intData[0] <<std::endl;
  std::cout << "Number of constraints: " << _intData[1] <<std::endl;
  std::cout << "itol, itask, istate, iopt, lrw, liw, jt: (for details on what are these variables see opkdmain.f)" <<std::endl;
  std::cout << _intData[2] << ", " << _intData[3] << ", " << _intData[4] << ", " << _intData[5] << ", " << _intData[6]  << ", " << _intData[7]  << ", " << _intData[8] <<std::endl;
  std::cout << "====================================" <<std::endl;
}
