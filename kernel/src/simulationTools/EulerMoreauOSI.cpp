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
#include "EulerMoreauOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderType2R.hpp"
#include "FirstOrderType1R.hpp"
#include "NonSmoothLaw.hpp"
#include "CxxStd.hpp"
#include "OneStepNSProblem.hpp"
#include "BlockVector.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace RELATION;

// --- constructor with theta parameter value  ---
EulerMoreauOSI::EulerMoreauOSI(double theta):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
  _steps=1;
  _theta = theta;
}

// --- constructor from a set of data ---
EulerMoreauOSI::EulerMoreauOSI(double theta, double gamma):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _useGammaForRelation(false)
{
  _levelMinForOutput= 0;
  _levelMaxForOutput =0;
  _levelMinForInput =0;
  _levelMaxForInput =0;
  _steps=1;
  _theta = theta;
  _gamma = gamma;
  _useGamma = true;
}

const SimpleMatrix EulerMoreauOSI::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "EulerMoreauOSI::getW(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "EulerMoreauOSI::getW(ds): W[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W); // Copy !!
}

SP::SimpleMatrix EulerMoreauOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::W(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
}


const SimpleMatrix EulerMoreauOSI::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "EulerMoreauOSI::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions
         && "EulerMoreauOSI::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions); // Copy !!
}

SP::SiconosMatrix EulerMoreauOSI::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::WBoundaryConditions(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions;
}

void EulerMoreauOSI::initializeDynamicalSystem(Model& m, double t, SP::DynamicalSystem ds)
{
  // Get work buffers from the graph
  const DynamicalSystemsGraph::VDescriptor& dsv = _dynamicalSystemsGraph->descriptor(ds);
  VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(dsv).workVectors;

  // W initialization
  initializeIterationMatrixW(t, ds, dsv);

  // buffers allocation (into the graph)
  workVectors.resize(OneStepIntegrator::work_vector_of_vector_size);
  workVectors[OneStepIntegrator::residu].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::residu_free].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::free].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::x_partial_ns].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::delta_x_for_relation].reset(new SiconosVector(ds->dimension()));
  workVectors[OneStepIntegrator::local_buffer].reset(new SiconosVector(ds->dimension()));
}

void EulerMoreauOSI::initializeInteraction(double t0, Interaction &inter,
					   InteractionProperties& interProp,
					   DynamicalSystemsGraph & DSG)
{
  SP::DynamicalSystem ds1= interProp.source;
  SP::DynamicalSystem ds2= interProp.target;
  assert(ds1);
  assert(ds2);

  assert(interProp.DSlink);

  VectorOfBlockVectors& DSlink = *interProp.DSlink;

  Relation &relation =  *inter.relation();
  RELATION::TYPES relationType = relation.getType();

  /* Check that the interaction has the correct initialization for y and lambda */

  bool isInitializationNeeded = false;
  if (!(inter.lowerLevelForOutput() <= _levelMinForOutput && inter.upperLevelForOutput()  >= _levelMaxForOutput ))
  {
    //RuntimeException::selfThrow("EulerMoreauOSI::initializeInteraction, we must resize _y");
    inter.setLowerLevelForOutput(_levelMinForOutput);
    inter.setUpperLevelForOutput(_levelMaxForOutput );
    isInitializationNeeded = true;
  }
  if (!(inter.lowerLevelForInput() <= _levelMinForInput && inter.upperLevelForInput() >= _levelMaxForInput ))
  {
    // RuntimeException::selfThrow("EulerMoreauOSI::initializeInteraction, we must resize _lambda");
    inter.setLowerLevelForInput(_levelMinForInput);
    inter.setUpperLevelForInput(_levelMaxForInput);
    isInitializationNeeded = true;
  }
  if (isInitializationNeeded)
    inter.init();

  bool computeResidu = relation.requireResidu();
  inter.initializeMemory(computeResidu,_steps);
  if(checkOSI(DSG.descriptor(ds1)))
    {
      DEBUG_PRINTF("ds1->number() %i is taken in to account\n", ds1->number());
      assert(DSG.properties(DSG.descriptor(ds1)).workVectors);
      VectorOfVectors &workVds1 = *DSG.properties(DSG.descriptor(ds1)).workVectors;


      if (relationType == FirstOrder)
	{
	  if (!DSlink[FirstOrderR::xfree])
	    {
	      DSlink[FirstOrderR::xfree].reset(new BlockVector());
	      DSlink[FirstOrderR::xfree]->insertPtr(workVds1[OneStepIntegrator::free]);
	    }
	  else
	    DSlink[FirstOrderR::xfree]->setVectorPtr(0,workVds1[OneStepIntegrator::free]);

	  if (!DSlink[FirstOrderR::xPartialNS])
	    {
	      DSlink[FirstOrderR::xPartialNS].reset(new BlockVector());
	      DSlink[FirstOrderR::xPartialNS]->insertPtr(workVds1[OneStepIntegrator::x_partial_ns]);
	    }
	  else
	    DSlink[FirstOrderR::xPartialNS]->setVectorPtr(0,workVds1[OneStepIntegrator::x_partial_ns]);
	  if (!DSlink[FirstOrderR::deltax])
	    {
	      DSlink[FirstOrderR::deltax].reset(new BlockVector());
	      DSlink[FirstOrderR::deltax]->insertPtr(workVds1[OneStepIntegrator::delta_x_for_relation]);
	    }
	  else
	    DSlink[FirstOrderR::deltax]->setVectorPtr(0,workVds1[OneStepIntegrator::delta_x_for_relation]);
	}
    }
  DEBUG_PRINTF("ds1->number() %i\n",ds1->number());
  DEBUG_PRINTF("ds2->number() %i\n",ds2->number());

  if (ds1 != ds2)
    {
      DEBUG_PRINT("ds1 != ds2\n");

      if(checkOSI(DSG.descriptor(ds2)))
	{
	  DEBUG_PRINTF("ds2->number() %i is taken in to account\n",ds2->number());
	  assert(DSG.properties(DSG.descriptor(ds2)).workVectors);
	  VectorOfVectors &workVds2 = *DSG.properties(DSG.descriptor(ds2)).workVectors;
	  if (relationType == FirstOrder)
	    {
	      if (!DSlink[FirstOrderR::xfree])
		{
		  DSlink[FirstOrderR::xfree].reset(new BlockVector());
		  //dummy insertion to reserve first vector for ds1
		  DSlink[FirstOrderR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
		  DSlink[FirstOrderR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);
		}
	      else
		DSlink[FirstOrderR::xfree]->insertPtr(workVds2[OneStepIntegrator::free]);

	      if (!DSlink[FirstOrderR::xPartialNS])
		{
		  DSlink[FirstOrderR::xPartialNS].reset(new BlockVector());
		  //dummy insertion to reserve first vector for ds1
		  DSlink[FirstOrderR::xPartialNS]->insertPtr(workVds2[OneStepIntegrator::x_partial_ns]);
		  DSlink[FirstOrderR::xPartialNS]->insertPtr(workVds2[OneStepIntegrator::x_partial_ns]);
		}
	      else
		DSlink[FirstOrderR::xPartialNS]->insertPtr(workVds2[OneStepIntegrator::x_partial_ns]);


	      if (!DSlink[FirstOrderR::deltax])
		{
		  DSlink[FirstOrderR::deltax].reset(new BlockVector());
		  //dummy insertion to reserve first vector for ds1
		  DSlink[FirstOrderR::deltax]->insertPtr(workVds2[OneStepIntegrator::delta_x_for_relation]);
		  DSlink[FirstOrderR::deltax]->insertPtr(workVds2[OneStepIntegrator::delta_x_for_relation]);
		}
	      else
		DSlink[FirstOrderR::deltax]->insertPtr(workVds2[OneStepIntegrator::delta_x_for_relation]);
	    }
	}
    }

  // Compute a first value for the output
  inter.computeOutput(t0, interProp, 0);

  // prepare the gradients
  relation.computeJach(t0, inter, interProp);
  for (unsigned int i = 0; i < inter.upperLevelForOutput() + 1; ++i)
    {
      inter.computeOutput(t0, interProp, i);
    }
  inter.swapInMemory();


}

void EulerMoreauOSI::initialize(Model& m)
{
  OneStepIntegrator::initialize(m);
  // Get initial time
  double t0 = _simulation->startingTime();

  std::cout << std::endl;
  // Compute W(t0) for all ds
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      initializeDynamicalSystem(m, t0,  ds);
    }

  SP::InteractionsGraph indexSet0 = m.nonSmoothDynamicalSystem()->topology()->indexSet0();
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
    {
      Interaction& inter = *indexSet0->bundle(*ui);

      initializeInteraction(t0, inter, indexSet0->properties(*ui), *_dynamicalSystemsGraph);
    }


}
void EulerMoreauOSI::initializeIterationMatrixW(double t, SP::DynamicalSystem ds, const DynamicalSystemsGraph::VDescriptor& dsv)
{
  // This function:
  // - allocate memory for a matrix W

  if(!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  if(_dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - W(ds) is already in the map and has been initialized.");


  unsigned int sizeW = ds->dimension(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W


  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);
  // 1 - First order non linear systems
  if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      //    // Memory allocation for W

      // W =  M - h*_theta* [jacobian_x f(t,x,z)]
      SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS> (ds);

      // Copy M or I if M is Null into W
      if(d->M())
	//      *W = *d->M();
	_dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->M()));

      else
	{
	  //W->eye();
	  _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(sizeW, sizeW)); // Warning if the Jacobian is a sparse matrix
	  _dynamicalSystemsGraph->properties(dsv).W->eye();
	}

      // d->computeJacobianfx(t); // Computation of JacxF is not required here
      // since it must have been done in OSI->initialize, before a call to this function.
      SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(dsv).W;
      // Add -h*_theta*jacobian_XF to W
      scal(-h * _theta, *d->jacobianfx(), *W, false);
    }
  else RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.




}


void EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if(!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions(t,ds) - ds == NULL");

  if(!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixW(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);


  RuntimeException::selfThrow("EulerMoreauOSI::initializeIterationMatrixWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeWBoundaryConditions(SP::DynamicalSystem ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initializeIterationMatrixWBoundaryConditions.

  assert(ds &&
         "EulerMoreauOSI::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);
  //unsigned int dsN = ds->number();
  RuntimeException::selfThrow("EulerMoreauOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeW(double t, DynamicalSystem& ds, DynamicalSystemsGraph::VDescriptor& dsv,
                              SiconosMatrix& W)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, W is supposed to exist and not to be null
  // Memory allocation has been done during initializeIterationMatrixW.

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(ds);

  // 1 - First order non linear systems
  if(dsType == Type::FirstOrderNonLinearDS)
    {
      // W =  M - h*_theta* [jacobian_x f(t,x,z)]
      FirstOrderNonLinearDS& d = static_cast<FirstOrderNonLinearDS&>(ds);

      // Copy M or I if M is Null into W
      if(d.M())
	W = *d.M();
      else
	W.eye();

      d.computeJacobianfx(t);
      // Add -h*_theta*jacobianfx to W
      scal(-h * _theta, *d.jacobianfx(), W, false);
    }
  // 2 - First order linear systems
  else if(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&>(ds);
      if(dsType == Type::FirstOrderLinearDS)
	d.computeA(t);

      if(d.M())
	W = *d.M();
      else
	W.eye();

      scal(-h * _theta, *d.A(), W, false);
    }
  else RuntimeException::selfThrow("EulerMoreauOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

  //  if (_useGamma)
  {

    InteractionsGraph& indexSet = *_simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0);

    DynamicalSystemsGraph::OEIterator oei, oeiend;
    InteractionsGraph::VDescriptor ivd;
    SP::SiconosMatrix K;
    SP::Interaction inter;
    for(std11::tie(oei, oeiend) = _dynamicalSystemsGraph->out_edges(dsv); oei != oeiend; ++oei)
      {
	inter = _dynamicalSystemsGraph->bundle(*oei);
	ivd = indexSet.descriptor(inter);
	FirstOrderR& rel = static_cast<FirstOrderR&>(*inter->relation());
	K = rel.K();
	if(!K) K = (*indexSet.properties(ivd).workMatrices)[FirstOrderR::mat_K];
	if(K)
	  {
	    scal(-h * _gamma, *K, W, false);
	  }
      }
  }
  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}



double EulerMoreauOSI::computeResidu()
{
  DEBUG_PRINT("EulerMoreauOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "EulerMoreauOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);


  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      ds = _dynamicalSystemsGraph->bundle(*dsi);
      dsType = Type::value(*ds); // Its type
      // XXX TMP hack -- xhub
      // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
      // Maurice will do that with subgraph :)
      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
      SiconosVector& residuFree = *workVectors[OneStepIntegrator::residu_free];
      SiconosVector& residu = *workVectors[OneStepIntegrator::residu];
      // 1 - First Order Non Linear Systems
      if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS)
	{
	  // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
	  //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
	  //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

	  // Note: indices k/k+1 corresponds to value at the beginning/end of the time step.
	  // Newton iterate are x and r

	  FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

	  // 1 - Compute the free residu (purely on the "smooth" dynamics)

	  residuFree = *d.x(); // last saved value for x: could be x_k or x_{k+1}^alpha
	  SiconosVector& xold = *d.xMemory()->getSiconosVector(0);
	  residuFree -= xold; // state x_k (at previous time step)

	  SP::SiconosMatrix M = d.M();
	  if(M)
	    prod(*M, residuFree, residuFree, true);
	  // at this step, we have residuFree = M(x - x_k)
	  DEBUG_PRINT("EulerMoreauOSI::computeResidu residuFree = M(x - x_k)\n");
	  DEBUG_EXPR(residuFree.display());

	  if(d.f())
	    {
	      double coef = -h * (1 - _theta);
	      if(dsType == Type::FirstOrderLinearDS)
		{
		  // computes f(t_k,x_k)
		  //This computation is done since fold not  is up to date.
		  d.computef(told, xold);
		  // residuFree += -h * (1 - _theta) * f(t_k,x_k)
		  scal(coef, *d.f(), residuFree, false);
		}
	      else // FirstOrderNonLinearDS
		{
		  // residuFree += -h * (1 - _theta) * f(t_k,x_k)
		  scal(coef, *d.fold(), residuFree, false);
		}

	      // computes f(t_{x+1}, x_{k+1}^alpha)
	      d.computef(t);
	      coef = -h * _theta;
	      // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
	      scal(coef, *d.f(), residuFree, false);
	    }

	  // now compute the residu = residuFree - h*gamma*r - h*(1-gamma)r_k
	  residu = residuFree;

	  if(!_useGamma)  // no gamma
	    {
	      scal(-h, *d.r(), residu, false); // residu = residu - h*r
	    }
	  else
	    {
	      scal(-h*_gamma, *d.r(), residu, false);
	      scal(-h*(1-_gamma), *d.rMemory()->getSiconosVector(0), residu, false);
	    }

	  normResidu = residu.norm2();
	  DEBUG_EXPR(residu.display());
	}
      // 2 - First Order Linear Systems with Time Invariant coefficients
      else if(dsType == Type::FirstOrderLinearTIDS)
	{
	  FirstOrderLinearTIDS& d = *std11::static_pointer_cast<FirstOrderLinearTIDS>(ds);
	  //Don't use W because it is LU factorized
	  //Residu : R_{free} = M(x^{\alpha}_{k+1} - x_{k}) -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1})
	  SP::SiconosVector b = d.b();
	  if(b)
	    residuFree = *b;
	  else
	    residuFree.zero();

	  SiconosVector& xBuffer = *workVectors[OneStepIntegrator::local_buffer];
	  SP::SiconosMatrix A = d.A();

	  if(A)  // residuFree += -h( A (\theta x_{k+1}^{\alpha} + (1-\theta) x_k)
	    {

	      prod(*A, *d.xMemory()->getSiconosVector(0), xBuffer, true);
	      double coef = -h * (1 - _theta);
	      scal(coef, xBuffer, residuFree, false);

	      prod(*A, *d.x(), xBuffer, true);
	      coef = -h * _theta;
	      scal(coef, xBuffer, residuFree, false);
	    }

	  // final touch, residuFree += M(x_{k+1}^{\alpha} - x_k)
	  xBuffer = *d.x() - *d.xMemory()->getSiconosVector(0);
	  SP::SiconosMatrix M = d.M();
	  if(M)
	    {
	      prod(*M, xBuffer, residuFree, false);
	    }
	  else
	    {
	      residuFree += xBuffer;
	    }


	}
      else
	RuntimeException::selfThrow("EulerMoreauOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

      DEBUG_PRINT("EulerMoreauOSI::computeResidu final residuFree\n");
      DEBUG_EXPR(residuFree.display());


      if(normResidu > maxResidu) maxResidu = normResidu;

    }
  return maxResidu;
}

void EulerMoreauOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_PRINT("EulerMoreauOSI::computeFreeState() starts\n");

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d.rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.



  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      ds = _dynamicalSystemsGraph->bundle(*dsi);

      // XXX TMP hack -- xhub
      // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
      // Maurice will do that with subgraph :)

      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      dsType = Type::value(*ds); // Its type
      SiconosMatrix& W = *_dynamicalSystemsGraph->properties(*dsi).W; // Its W EulerMoreauOSI matrix of iteration.

      // 1 - First Order Non Linear Systems
      if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
	{
	  // xfree =  x - W^{-1} (ResiduFree - h(1-gamma)*rold)
	  // with ResiduFree = = M(x - x_k) - h*theta*f(t_{k+1}, x) - h*(1-theta)*f(t_k, x_k)

	  // to be updated at current time: W, f
	  // fold is f at t_k
	  // not time dependant: M
	  FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

	  // Get state i (previous time step) from Memories -> var. indexed with "Old"
	  //    SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

	  SiconosVector& x = *d.x(); // x = x_k or x = x_{k+1}^{\alpha}
	  // xfree gets ResiduFree at first
	  SiconosVector& xfree = *workVectors[OneStepIntegrator::free];
	  xfree = *workVectors[OneStepIntegrator::residu_free];

	  DEBUG_PRINT("EulerMoreauOSI::computeFreeState xfree <- residuFree\n");
	  DEBUG_EXPR(xfree.display());

	  if(_useGamma)
	    {
	      SiconosVector& rold = *d.rMemory()->getSiconosVector(0);
	      double coeff = -h * (1 - _gamma);
	      scal(coeff, rold, xfree, false); //  xfree += -h(1-gamma)*rold
	    }


	  // At this point xfree = (ResiduFree - h(1-gamma)*rold)
	  // -> Solve WX = xfree and set xfree = X
	  W.PLUForwardBackwardInPlace(xfree);

	  // at this point, xfree = W^{-1} (ResiduFree - h(1-gamma)*rold)
	  // -> compute real xfree = x - W^{-1} (ResiduFree - h(1-gamma)*rold)
	  xfree *= -1.0;
	  xfree += x;

	  DEBUG_EXPR(xfree.display());

	  // now the crazy intermediate variables
	  // xPartialNS was updated before this fonction call
	  // It constains either 0 (first Newton iterate)
	  // or g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
	  SiconosVector& xPartialNS = *workVectors[OneStepIntegrator::x_partial_ns];
	  DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS from Interaction\n");
	  DEBUG_EXPR(xPartialNS.display());

	  // -> Solve WX = g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
	  // and set xPartialNS = X
	  W.PLUForwardBackwardInPlace(xPartialNS);
	  scal(h, xPartialNS, xPartialNS);

	  // compute real xPartialNS = xfree + ...
	  xPartialNS += xfree;
	  DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS real value\n");
	  DEBUG_EXPR(xPartialNS.display());

	  // deltaxForRelation = (\widetilde{K}_{k+1}^{\alpha})^{-1} xPartialNS - x
	  SiconosVector& deltaxForRelation = *workVectors[OneStepIntegrator::delta_x_for_relation];
	  deltaxForRelation = xPartialNS;

	  deltaxForRelation -= x;

	  DEBUG_EXPR(deltaxForRelation.display());

	  // have a look at the end of the DevNotes for this part
	  if(_useGammaForRelation)
	    {
	      if(!(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS))
		RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - _useGammaForRelation == true is only implemented for FirstOrderLinearDS or FirstOrderLinearTIDS");

	      deltaxForRelation = xfree;

	      scal(_gamma, deltaxForRelation, deltaxForRelation);
	      SiconosVector& xold = *d.xMemory()->getSiconosVector(0);

	      scal(1.0 - _gamma, xold, deltaxForRelation, false);
	    }

	  // some output
	  DEBUG_EXPR(xfree.display(););
	  DEBUG_EXPR(xPartialNS.display(););
	  DEBUG_EXPR(deltaxForRelation.display(););

	}
      else
	RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
}

void EulerMoreauOSI::prepareNewtonIteration(double time)
{
  // XXX TMP hack -- xhub
  // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
  // Maurice will do that with subgraph :)
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      DynamicalSystemsGraph::VDescriptor dsv = _dynamicalSystemsGraph->descriptor(ds);
      SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(*dsi).W;
      computeW(time, *ds, dsv, *W);
    }
}

/// @cond

struct EulerMoreauOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) :
    _osnsp(p), _inter(inter) {};

  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};

/// @endcond

void EulerMoreauOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */

  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  VectorOfVectors& workV = *indexSet->properties(vertex_inter).workVectors;
  VectorOfSMatrices& workM = *indexSet->properties(vertex_inter).workMatrices;
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();

  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;



  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix  C;
  SP::SiconosMatrix  D;
  SP::SiconosMatrix  F;
  SP::BlockVector deltax;
  SP::BlockVector Xfree;

  SP::SiconosVector H_alpha;

  deltax = DSlink[FirstOrderR::deltax];
  SiconosVector& yForNSsolver = *inter->yForNSsolver();

  Xfree = DSlink[FirstOrderR::xfree];

  assert(Xfree);


  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if(relationType == FirstOrder && (relationSubType == Type2R || relationSubType == NonLinearR))
    {
      SiconosVector& lambda = *inter->lambda(0);
      FirstOrderR& rel = *std11::static_pointer_cast<FirstOrderR>(mainInteraction->relation());
      C = rel.C();
      if(!C) C = workM[FirstOrderR::mat_C];
      D = rel.D();
      if(!D) D = workM[FirstOrderR::mat_D];

      if(D)
	{
	  coord[3] = D->size(1);
	  coord[5] = D->size(1);
	  subprod(*D, lambda, yForNSsolver, coord, true);

	  yForNSsolver *= -1.0;
	}
      else
	{
	  subscal(0, yForNSsolver, yForNSsolver, coord, true);
	}

      if(C)
	{
	  coord[3] = C->size(1);
	  coord[5] = C->size(1);
	  subprod(*C, *deltax, yForNSsolver, coord, false);

	}

      if(_useGammaForRelation)
	{
	  RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Type2R and H_alpha->getValue() should return the mid-point value");
	}

      H_alpha = inter->Halpha();
      assert(H_alpha);
      yForNSsolver += *H_alpha;
    }
  else if(relationType == FirstOrder && relationSubType == Type1R)
    {
      FirstOrderType1R& rel = *std11::static_pointer_cast<FirstOrderType1R>(mainInteraction->relation());
      C = rel.C();
      if(!C) C = workM[FirstOrderR::mat_C];
      F = rel.F();
      if(!F) F = workM[FirstOrderR::mat_F];
      assert(Xfree);
      assert(deltax);

      if(F)
	{
	  coord[3] = F->size(1);
	  coord[5] = F->size(1);
	  subprod(*F, *DSlink[FirstOrderR::z], yForNSsolver, coord, true);

	}
      if(C)
	{
	  coord[3] = C->size(1);
	  coord[5] = C->size(1);
	  subprod(*C, *Xfree, yForNSsolver, coord, false);

	}

      if(_useGammaForRelation)
	{
	  RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
	}
      H_alpha = inter->Halpha();
      if(H_alpha)
	{
	  yForNSsolver += *H_alpha;
	}
    }
  else // First Order Linear Relation
    {
      C = mainInteraction->relation()->C();
      if(!C) C = workM[FirstOrderR::mat_C];

      if(C)
	{

	  assert(Xfree);
	  assert(deltax);

	  coord[3] = C->size(1);
	  coord[5] = C->size(1);

	  if(_useGammaForRelation)
	    {
	      subprod(*C, *deltax, yForNSsolver, coord, true);
	    }
	  else
	    {
	      subprod(*C, *Xfree, yForNSsolver, coord, true);
	    }
	}

      if(relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
	{
	  // In the first order linear case it may be required to add e + FZ to y.
	  // y = CXfree + e + FZ
	  SP::SiconosVector e;
	  if(relationSubType == LinearTIR)
	    {
	      e = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->e();
	      F = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->F();
	    }
	  else
	    {
	      e = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->e();
	      if(!e) e = workV[FirstOrderR::e];
	      F = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->F();
	      if(!F) F = workM[FirstOrderR::mat_F];
	    }

	  if(e)
	    yForNSsolver += *e;

	  if(F)
	    {
	      coord[3] = F->size(1);
	      coord[5] = F->size(1);
	      subprod(*F, *DSlink[FirstOrderR::z], yForNSsolver, coord, false);
	    }
	}

    }
}

void EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  //double h = tend - tinit;
  tout = tend;


  SP::SiconosMatrix W;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      Type::Siconos dsType = Type::value(*ds);
      RuntimeException::selfThrow("EulerMoreauOSI::integrate - not yet implemented for Dynamical system type :" + dsType);
    }
}

void EulerMoreauOSI::updateState(const unsigned int level)
{

  DEBUG_PRINT("EulerMoreauOSI::updateState\n");

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if(useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;

  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

      // Get the DS type
      Type::Siconos dsType = Type::value(*ds);
      VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(*dsi).workVectors;

      SimpleMatrix&  W = *_dynamicalSystemsGraph->properties(*dsi).W;

      if(dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
	{
	  FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
	  SiconosVector& x = *ds->x();
	  DEBUG_PRINT("EulerMoreauOSI::updateState Old value of x\n");
	  DEBUG_EXPR(x.display());
	  DEBUG_PRINT("EulerMoreauOSI::updateState residu value\n");
	  DEBUG_EXPR(d.r()->display());

	  // TODO ???
	  bool baux = (useRCC && dsType == Type::FirstOrderNonLinearDS && _simulation->relativeConvergenceCriterionHeld());
	  if(level != LEVELMAX)
	    {

	      //    SP::SiconosVector xFree = d->xFree();

	      // Save value of q in local_buffer for relative convergence computation
	      if(baux)
		*workVectors[OneStepIntegrator::local_buffer] = x;

	      //        std::cout <<boolalpha << _useGamma << std::endl;
	      //        std::cout <<_gamma << std::endl;
	      if(_useGamma)
		{
		  // XXX UseGamma broken ? -- xhub
		  scal(_gamma * h, *d.r(), x); // x = gamma*h*r
		}
	      else
		{
		  scal(h, *d.r(), x); // x = h*r
		}

	      W.PLUForwardBackwardInPlace(x); // x = h* W^{-1} *r

	      x += *workVectors[OneStepIntegrator::free]; // x+=xfree
	    }
	  else
	    {
	      RuntimeException::selfThrow("EulerMoreauOSI::updateState - level != LEVELMAX is not supposed to happen !");
	      x = *workVectors[OneStepIntegrator::free]; // x = xfree
	    }

	  if(baux)
	    {
	      *workVectors[OneStepIntegrator::local_buffer] -= x;
	      double aux = (workVectors[OneStepIntegrator::local_buffer]->norm2()) / (ds->normRef());
	      if(aux > RelativeTol)
		_simulation->setRelativeConvergenceCriterionHeld(false);
	    }
	  DEBUG_PRINT("EulerMoreauOSI::updateState New value of x\n");
	  DEBUG_EXPR(x.display());
	}
      else RuntimeException::selfThrow("EulerMoreauOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
    }
}

void EulerMoreauOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======" <<std::endl;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for(std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
    {
      if(!checkOSI(dsi)) continue;
      SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
      std::cout << "--------------------------------" <<std::endl;
      std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
      if(_dynamicalSystemsGraph->properties(*dsi).W) _dynamicalSystemsGraph->properties(*dsi).W->display();
      else std::cout << "-> NULL" <<std::endl;
      std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
    }
  std::cout << "================================" <<std::endl;
}
