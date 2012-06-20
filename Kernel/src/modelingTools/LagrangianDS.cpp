
/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "LagrangianDS.hpp"
#include "LagrangianDSXML.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"
#include <iostream>
using namespace std;

// Private function to set linked with members of Dynamical top class
void LagrangianDS::connectToDS()
{
  // dim
  _n = 2 * _ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianfx.

  // Initial conditions
  _x0.reset(new SiconosVector(*_q0, *_velocity0));

  // Current State: \f$ x \f$ and rhs = \f$ \dot x \f$

  _x[0].reset(new SiconosVector(*_q[0], *_q[1]));
  _x[1].reset(new SiconosVector(*_q[1], *_q[2]));
  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
}
void LagrangianDS::zeroPlugin()
{
  DynamicalSystem::zeroPlugin();
  _pluginMass.reset(new PluggedObject());
  _pluginFInt.reset(new PluggedObject());
  _pluginFExt.reset(new PluggedObject());
  _pluginNNL.reset(new PluggedObject());
  _pluginJacqFInt.reset(new PluggedObject());
  _pluginJacqDotFInt.reset(new PluggedObject());
  _pluginJacqNNL.reset(new PluggedObject());
  _pluginJacqDotNNL.reset(new PluggedObject());
}

LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0):
  DynamicalSystem(2 * newQ0->size()), _ndof(newQ0->size())
{
  zeroPlugin();
  // -- Memory allocation for vector and matrix members --
  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SiconosVector(*_q0));
  _q[1].reset(new SiconosVector(*_velocity0));
  _q[2].reset(new SiconosVector(_ndof));
  _residuFree.reset(new SiconosVector(getDim()));
  //   _xp.reset(new SiconosVector(getDim()));
  //   _xq.reset(new SiconosVector(getDim()));
  //   mXfree.reset(new SiconosVector(getDim()));
  //   r.reset(new SiconosVector(getDim()));

  // set allocation flags: true for required input, false for others
  _p.resize(3);
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(Type::LagrangianDS), _ndof(0)
{
  zeroPlugin();
  // Protected constructor - Only call from derived class(es).
  _q.resize(3);
  _p.resize(3);
  // !!! No plug-in connection !!!
}

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(SP::DynamicalSystemXML dsxml):
  DynamicalSystem(dsxml), _ndof(0)
{
  zeroPlugin();
  // -- Lagrangian  xml object --
  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);

  // === Initial conditions ===
  // Warning: ndof is given by q0.size() !!
  if (! lgptr->hasQ0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, q0 is a required input");

  _q0.reset(new SiconosVector(lgptr->getQ0())); // required node in xml file
  _ndof = _q0->size();

  if (! lgptr->hasVelocity0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, v0 is a required input");

  _velocity0.reset(new SiconosVector(lgptr->getVelocity0())); // required node in xml file

  if (_velocity0->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS::xml constructor - size of input velocity0 differs from ndof");

  // --- Current State (optional input) ---

  _q.resize(3);

  if (lgptr->hasQ())
    _q[0].reset(new SiconosVector(lgptr->getQ())); // Means that first q is different from q0 ?? Strange case ...
  else
    _q[0].reset(new SiconosVector(*_q0));           // initialize q with q0
  if (lgptr->hasVelocity())
    _q[1].reset(new SiconosVector(lgptr->getVelocity0())); // same remark as for q
  else
    _q[1].reset(new SiconosVector(*_velocity0));

  _q[2].reset(new SiconosVector(_ndof));
  _residuFree.reset(new SiconosVector(getDim()));
  _p.resize(3);
  // Memories
  if (lgptr->hasQMemory())   // qMemory
    _qMemory.reset(new SiconosMemory(lgptr->getQMemoryXML()));

  if (lgptr->hasVelocityMemory())   // velocityMemory
    _velocityMemory.reset(new SiconosMemory(lgptr->getVelocityMemoryXML()));

  string plugin;
  // mass
  if (lgptr->isMassPlugin()) // if mass is plugged
  {
    plugin = lgptr->getMassPlugin();
    setComputeMassFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
  }
  else
    _mass.reset(new SimpleMatrix(lgptr->getMassMatrix()));

  // === Optional inputs ===

  // fInt
  if (lgptr->hasFInt())  // if fInt is given
  {
    if (lgptr->isFIntPlugin()) // if fInt is plugged
    {
      plugin = lgptr->getFIntPlugin();
      setComputeFIntFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _fInt.reset(new SiconosVector(lgptr->getFIntVector()));
  }

  // _fExt
  if (lgptr->hasFExt())  // if _fExt is given
  {
    if (lgptr->isFExtPlugin())// if _fExt is plugged
    {
      plugin = lgptr->getFExtPlugin();
      setComputeFExtFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _fExt.reset(new SiconosVector(lgptr->getFExtVector()));
  }

  // NNL
  if (lgptr ->hasNNL())// if NNL is given
  {
    if (lgptr->isNNLPlugin())// if NNL is plugged
    {
      plugin = lgptr->getNNLPlugin();
      setComputeNNLFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _NNL.reset(new SiconosVector(lgptr->getNNLVector()));
  }

  // jacobian q of fInt
  if (lgptr ->hasJacobianFInt(0))// if given
  {
    if (lgptr->isJacobianFIntPlugin(0))// if is plugged
    {
      plugin = lgptr->getJacobianFIntPlugin(0);
      setComputeJacobianFIntqFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _jacobianFIntq.reset(new SimpleMatrix(lgptr->getJacobianFIntMatrix(0)));
  }

  // jacobian qdot of fInt
  if (lgptr ->hasJacobianFInt(1))// if given
  {
    if (lgptr->isJacobianFIntPlugin(1))// if is plugged
    {
      plugin = lgptr->getJacobianFIntPlugin(1);
      setComputeJacobianFIntqDotFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _jacobianFIntqDot.reset(new SimpleMatrix(lgptr->getJacobianFIntMatrix(1)));
  }


  // Jacobian q of NNL
  if (lgptr -> hasJacobianNNL(0)) // if given
  {
    if (lgptr->isJacobianNNLPlugin(0))// if is plugged
    {
      plugin = lgptr->getJacobianNNLPlugin(0);
      setComputeJacobianNNLqFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _jacobianNNLq.reset(new SimpleMatrix(lgptr->getJacobianNNLMatrix(0)));
  }
  // Jacobian qdot of NNL
  if (lgptr -> hasJacobianNNL(1)) // if given
  {
    if (lgptr->isJacobianNNLPlugin(1))// if is plugged
    {
      plugin = lgptr->getJacobianNNLPlugin(1);
      setComputeJacobianNNLqDotFunction(SSL::getPluginName(plugin), SSL::getPluginFunctionName(plugin));
    }
    else
      _jacobianNNLqDot.reset(new SimpleMatrix(lgptr->getJacobianNNLMatrix(1)));
  }

}

// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum Lagrangian System form: \f$ M\ddot q = p \f$
/**/
LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0, SP::SiconosMatrix newMass):
  DynamicalSystem(2 * newQ0->size()), _ndof(newQ0->size())
{
  zeroPlugin();
  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  // Mass matrix
  _mass = newMass;

  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SiconosVector(*_q0));
  _q[1].reset(new SiconosVector(*_velocity0));
  _q[2].reset(new SiconosVector(_ndof));
  _residuFree.reset(new SiconosVector(getDim()));


  _p.resize(3);
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0, const string& massName):
  DynamicalSystem(), _ndof(newQ0->size())
{
  zeroPlugin();
  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SiconosVector(*_q0));
  _q[1].reset(new SiconosVector(*_velocity0));
  _q[2].reset(new SiconosVector(_ndof));
  _residuFree.reset(new SiconosVector(getDim()));

  // Mass
  setComputeMassFunction(SSL::getPluginName(massName), SSL::getPluginFunctionName(massName));

  _p.resize(3);
}

// Destructor
LagrangianDS::~LagrangianDS()
{
}

bool LagrangianDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (_ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0
  if (! _q0 || ! _velocity0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (! _mass)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - Mass not set.");
    output = false;
  }

  // fInt
  if ((_fInt && _pluginFInt->fPtr) && (! _jacobianFIntq || ! _jacobianFIntqDot))
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
    output = false;
  }

  // NNL
  if ((_NNL  && _pluginNNL->fPtr) && (! _jacobianNNLq || ! _jacobianNNLqDot))
    // ie if NNL is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined NNL but not its Jacobian (according to q and velocity).");
    output = false;
  }

  if (!output) cout << "LagrangianDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void LagrangianDS::initializeNonSmoothInput(unsigned int level)
{
  if (!_p[level])
    _p[level].reset(new SiconosVector(_ndof));
}

void LagrangianDS::initForces()
{

  _forces.reset(new SiconosVector(_ndof));

  _jacobianqForces.reset(new SimpleMatrix(_ndof, _ndof));
  _jacobianqDotForces.reset(new SimpleMatrix(_ndof, _ndof));
}

void LagrangianDS::initRhs(double time)
{
  _workMatrix.resize(sizeWorkMat);

  // Solve Mq[2]=fL+p.
  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_forces)
  {
    computeForces(time);
    *_q[2] += *_forces;
  }
  computeMass();
  // Copy of Mass into _workMatrix for LU-factorization.

  _workMatrix[invMass].reset(new SimpleMatrix(*_mass));
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);

  bool flag1 = false, flag2 = false;
  if (_jacobianqForces)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqForces(time);

    _workMatrix[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianqForces));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc10]);
    flag1 = true;
  }

  if (_jacobianqDotForces)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianqDotForces(time);
    _workMatrix[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianqDotForces));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc11]);
    flag2 = true;
  }

  _workMatrix[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::ZERO));
  _workMatrix[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, Siconos::IDENTITY));

  if (flag1 && flag2)
    _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[jacobianXBloc11]));
  else if (flag1) // flag2 = false
    _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[zeroMatrix]));
  else if (flag2) // flag1 = false
    _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[jacobianXBloc11]));
  else
    _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[zeroMatrix]));
}

void LagrangianDS::initialize(double time, unsigned int sizeOfMemory)
{

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *_q[0] = *_q0;
  *_q[1] = *_velocity0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SiconosVector(1));
  if (_pluginNNL->fPtr && !_NNL)
    _NNL.reset(new SiconosVector(_ndof));
  if (_pluginJacqDotNNL->fPtr && !_jacobianNNLqDot)
    _jacobianNNLqDot.reset(new SimpleMatrix(_ndof, _ndof));
  if (_pluginJacqNNL->fPtr && ! _jacobianNNLq)
    _jacobianNNLq.reset(new SimpleMatrix(_ndof, _ndof));

  if (_pluginFExt->fPtr && !_fExt)
    _fExt.reset(new SiconosVector(_ndof));

  if (_pluginFInt->fPtr && ! _fInt)
    _fInt.reset(new SiconosVector(_ndof));
  if (_pluginJacqFInt->fPtr && !_jacobianFIntq)
    _jacobianFIntq.reset(new SimpleMatrix(_ndof, _ndof));
  if (_pluginJacqDotFInt->fPtr && !_jacobianFIntqDot)
    _jacobianFIntqDot.reset(new SimpleMatrix(_ndof, _ndof));





  //
  if (!_workFree)
    _workFree.reset(new SiconosVector(getDim()));
  // Memory allocation for fL and its jacobians.
  initForces();


  if (_boundaryConditions)
  {
    _reactionToBoundaryConditions.reset(new SiconosVector(_boundaryConditions->velocityIndices()->size()));
  }


  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  if (!_mass)
    _mass.reset(new SimpleMatrix(_ndof, _ndof));



  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  //initRhs(time);




}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if (newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if (! _q[0])
    _q[0].reset(new SiconosVector(newValue));
  else
    *_q[0] = newValue;
}

void LagrangianDS::setQPtr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");
  _q[0] = newPtr;

}

void LagrangianDS::setQ0(const SiconosVector& newValue)
{
  if (newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if (! _q0)
    _q0.reset(new SiconosVector(newValue));
  else
    *_q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");
  _q0 = newPtr;
}

void LagrangianDS::setVelocity(const SiconosVector& newValue)
{
  if (newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (! _q[1])
    _q[1].reset(new SiconosVector(newValue));
  else
    *_q[1] = newValue;
}

void LagrangianDS::setVelocityPtr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");
  _q[1] = newPtr;
}


void LagrangianDS::setVelocity0Ptr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");
  _velocity0 = newPtr;
}

SP::SiconosVector LagrangianDS::acceleration() const
{
  return _q[2];
}

void LagrangianDS::setP(const SiconosVector& newValue, unsigned int level)
{
  if (newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if (! _p[level])
    _p[level].reset(new SiconosVector(newValue));
  else
    *(_p[level]) = newValue;
}

void LagrangianDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
{

  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");
  _p[level] = newPtr;
}

void LagrangianDS::computeMass()
{
  if (_pluginMass->fPtr)
    ((FPtrMass)_pluginMass->fPtr)(_ndof, &(*_q[0])(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeMass(SP::SiconosVector q2)
{
  if (_pluginMass->fPtr)
    ((FPtrMass)_pluginMass->fPtr)(_ndof, &(*q2)(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFInt(double time)
{
  if (_pluginFInt->fPtr)
    ((FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginFInt->fPtr)
    ((FPtr6)_pluginFInt->fPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFExt(double time)
{
  if (_pluginFExt->fPtr)
    ((FPtrFExt)_pluginFExt->fPtr)(time, _ndof, &(*_fExt)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeNNL()
{
  if (_pluginNNL->fPtr)
    ((FPtr5)_pluginNNL->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginNNL->fPtr)
    ((FPtr5)_pluginNNL->fPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianFIntq(double time)
{
  if (_pluginJacqFInt->fPtr)
    ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFIntq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFIntqDot(double time)
{
  if (_pluginJacqDotFInt->fPtr)
    ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFIntqDot)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZFInt(double time){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianFIntq(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginJacqFInt->fPtr)
    ((FPtr6)_pluginJacqFInt->fPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianFIntq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianFIntqDot(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginJacqDotFInt->fPtr)
    ((FPtr6)_pluginJacqDotFInt->fPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianFIntqDot)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianNNLq()
{
  if (_pluginJacqNNL->fPtr)
    ((FPtr5)_pluginJacqNNL->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianNNLq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianNNLqDot()
{
  if (_pluginJacqDotNNL->fPtr)
    ((FPtr5)_pluginJacqDotNNL->fPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianNNLqDot)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZNNL(){
//   if(computeJacobianZNNLPtr)
//       (computeJacobianZNNLPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianNNLq(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginJacqNNL->fPtr)
    ((FPtr5)_pluginJacqNNL->fPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianNNLq)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianNNLqDot(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (_pluginJacqDotNNL->fPtr)
    ((FPtr5)_pluginJacqDotNNL->fPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianNNLqDot)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZNNL(unsigned int i, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZNNLPtr)
//     (computeJacobianZNNLPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_forces)
  {
    computeForces(time);
    *_q[2] += *_forces;
  }

  // mass and inv(mass) computatiton
  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
    computeMass();


  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into _workMatrix[invMass].
  // -- Case 2: mass is not constant, we copy it into _workMatrix[invMass]
  // Then we proceed with PLUForwardBackward.

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  *_workMatrix[invMass] = *_mass;

  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);

}

void LagrangianDS::computeJacobianRhsx(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  if (!isDSup)
    computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children
  *_workMatrix[invMass] = *_mass;

  if (_jacobianqForces)
  {
    SP::SiconosMatrix bloc10 = _jacxRhs->block(1, 0);
    computeJacobianqForces(time);
    *bloc10 = *_jacobianqForces;
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc10);
  }

  if (_jacobianqDotForces)
  {
    SP::SiconosMatrix bloc11 = _jacxRhs->block(1, 1);
    computeJacobianqDotForces(time);
    *bloc11 = *_jacobianqDotForces;
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeForces(double time)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_forces)
  {
    // 1 - Computes the required functions
    computeFInt(time);
    computeFExt(time);
    computeNNL();

    // 2 - set fL = fExt - fInt - NNL

    // seems ok.
    if (_forces.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      _forces->zero();

      if (_fInt)
        *_forces -= *_fInt;

      if (_fExt)
        *_forces += *_fExt;

      if (_NNL)
        *_forces -= *_NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_forces)
  {
    // 1 - Computes the required functions
    computeFInt(time, q2, v2);
    computeFExt(time);
    computeNNL(q2, v2);

    // seems ok.
    if (_forces.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      _forces->zero();

      if (_fInt)
        *_forces -= *_fInt;

      if (_fExt)
        *_forces += *_fExt;

      if (_NNL)
        *_forces -= *_NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeJacobianqForces(double time)
{
  if (_jacobianqForces)
  {
    computeJacobianFIntq(time);
    computeJacobianNNLq();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianqForces->zero();
      if (_jacobianFIntq)
        *_jacobianqForces -= *_jacobianFIntq;
      if (_jacobianNNLq)
        *_jacobianqForces -= *_jacobianNNLq;
    }
  }
  //else nothing.
}
void LagrangianDS::computeJacobianqDotForces(double time)
{
  if (_jacobianqDotForces)
  {
    computeJacobianFIntqDot(time);
    computeJacobianNNLqDot();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianqDotForces->zero();
      if (_jacobianFIntqDot)
        *_jacobianqDotForces -= *_jacobianFIntqDot;
      if (_jacobianNNLqDot)
        *_jacobianqDotForces -= *_jacobianNNLqDot;
    }
  }
  //else nothing.
}
// void LagrangianDS::computeJacobianZFL( double time){
//    RuntimeException::selfThrow("LagrangianDS::computeJacobianZFL - not implemented");
// }

void LagrangianDS::saveSpecificDataToXML()
{
  // --- other data ---
  if (!_dsxml)
    RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DynamicalSystemXML does not exist");

  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(_dsxml);
  lgptr->setMassPlugin(_pluginMass->getPluginName());
  lgptr->setQ(*_q[0]);
  lgptr->setQ0(*_q0);
  lgptr->setQMemory(*_qMemory);
  lgptr->setVelocity(*_q[1]);
  lgptr->setVelocity0(*_velocity0);
  lgptr->setVelocityMemory(*_velocityMemory);

  // FExt
  if (lgptr->hasFExt())
  {
    if (!lgptr->isFExtPlugin())
      lgptr->setFExtVector(*_fExt);
  }
  else
    lgptr->setFExtPlugin(_pluginFExt->getPluginName());

  // FInt
  if (lgptr->hasFInt())
  {
    if (!lgptr->isFIntPlugin())
    {
      if (_fInt->size() > 0)
        lgptr->setFIntVector(*_fInt);
      else cout << "Warning : FInt can't be saved, the FInt vector is not defined." << endl;
    }
  }
  else
    lgptr->setFIntPlugin(_pluginFInt->getPluginName());

  for (unsigned int i = 0; i < 2; ++i)
  {
    // Jacobian FInt
    if (lgptr->hasJacobianFInt(i))
    {
      if (!lgptr->isJacobianFIntPlugin(i))
        lgptr->setJacobianFIntMatrix(i, i ? *_jacobianFIntqDot : *_jacobianFIntq);
    }
    else
      lgptr->setJacobianFIntPlugin(i, i ? _pluginJacqDotFInt->getPluginName() : _pluginJacqFInt->getPluginName());

    // JacobianNNLq
    if (lgptr->hasJacobianNNL(i))
    {
      if (!lgptr->isJacobianNNLPlugin(i))
        lgptr->setJacobianNNLMatrix(i, i ? *_jacobianNNLqDot : *_jacobianNNLq);
    }
    else
      lgptr->setJacobianNNLPlugin(i, i ? _pluginJacqDotNNL->getPluginName() : _pluginJacqNNL->getPluginName());
  }
  // NNL
  if (lgptr->hasNNL())
  {
    if (!lgptr->isNNLPlugin())
      lgptr->setNNLVector(*_NNL);
  }
  else
    lgptr->setNNLPlugin(_pluginNNL->getPluginName());
}

void LagrangianDS::display() const
{
  cout << "=====> Lagrangian System display (number: " << _number << ")." << endl;
  cout << "- _ndof : " << _ndof << endl;
  cout << "- q " << endl;
  if (_q[0]) _q[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- q0 " << endl;
  if (_q0) _q0->display();
  cout << "- v " << endl;
  if (_q[1]) _q[1]->display();
  else cout << "-> NULL" << endl;
  cout << "- v0 " << endl;
  if (_velocity0) _velocity0->display();
  else cout << "-> NULL" << endl;
  cout << "- p[0] " << endl;
  if (_p[0]) _p[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- p[1] " << endl;
  if (_p[1]) _p[1]->display();
  else cout << "-> NULL" << endl;
  cout << "- p[2] " << endl;
  if (_p[2]) _p[2]->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);
  if (steps == 0)
    cout << "Warning : LagragianDS::initMemory with size equal to zero" << endl;
  else
  {
    _qMemory.reset(new SiconosMemory(steps));
    _velocityMemory.reset(new SiconosMemory(steps));
    //swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{
  _xMemory->swap(_x[0]);
  _qMemory->swap(_q[0]);
  _velocityMemory->swap(_q[1]);
  // initialization of the reaction force due to the non smooth law
  // _p[1]->zero();
}

LagrangianDS* LagrangianDS::convert(DynamicalSystem* ds)
{
  LagrangianDS* lnlds = dynamic_cast<LagrangianDS*>(ds);
  return lnlds;
}
/*must be remove, replace by the RelativeConvergenceCriteron of the simulation*/

/*double LagrangianDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  SP::SiconosVector valRef = workV[NewtonSave];

  sub(*(q[0]),*valRef,*valRef);
  dsCvgIndic= valRef->norm2()/(valRef->norm2()+1);
  return (dsCvgIndic);
  }*/

// void LagrangianDS::computeqFree(double time, unsigned int level, SP::SiconosVector qFreeOut)
// {
//   // to compute qFree, derivative number level. Mainly used in EventDriven to compute second derivative
//   // of q for Output y computation.

//   if(qFreeOut->size()!=_ndof)
//     RuntimeException::selfThrow("LagrangianDS::computeqFree - Wrong size for output (different from _ndof)");


//   if(level!=2)
//     RuntimeException::selfThrow("LagrangianDS::computeqFree - Only implemented for second derivative.");

//   // Warning: we suppose that all other operators are up to date (FInt, FExt ...)

//   qFreeOut->zero();
//   if( _fInt )
//     *qFreeOut -= *_fInt;
//   if( _fExt )
//     *qFreeOut += *_fExt;
//   if( _NNL )
//     *qFreeOut -= *_NNL;

//   _workMatrix[invMass]->PLUForwardBackwardInPlace(*qFreeOut);
// }

void LagrangianDS::resetNonSmoothPart()
{
  if (_p[0])
    _p[0]->zero();
  if (_p[1])
    _p[1]->zero();
  if (_p[2])
    _p[2]->zero();
}

void LagrangianDS::resetNonSmoothPart(unsigned int level)
{
  if (_p[level])
    _p[level]->zero();
}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  SiconosVector tmp(*_p[1]);
  _workMatrix[invMass]->PLUForwardBackwardInPlace(tmp);
  *_q[1] += tmp;  // v+ = v- + p
}

void LagrangianDS::setComputeNNLFunction(const std::string& pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeNNLPtr, pluginPath,functionName);
  _pluginNNL->setComputeFunction(pluginPath, functionName);
}
void LagrangianDS::setComputeNNLFunction(FPtr5 fct)
{
  _pluginNNL->setComputeFunction((void *)fct);
  //    computeNNLPtr=fct;
}
void LagrangianDS::setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntqPtr, pluginPath,functionName);
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianFIntqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntqDotPtr, pluginPath,functionName);
  _pluginJacqDotFInt->setComputeFunction(pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianFIntqFunction(FPtr6 fct)
{
  _pluginJacqFInt->setComputeFunction((void *)fct);
}
void LagrangianDS::setComputeJacobianFIntqDotFunction(FPtr6 fct)
{
  _pluginJacqDotFInt->setComputeFunction((void *)fct);
}
void LagrangianDS::setComputeJacobianNNLqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqNNL->setComputeFunction(pluginPath, functionName); // Plugin::setFunction(&computeJacobianNNLqPtr, pluginPath,functionName);
}
void LagrangianDS::setComputeJacobianNNLqDotFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  _pluginJacqDotNNL->setComputeFunction(pluginPath, functionName);
  // Plugin::setFunction(&computeJacobianNNLqDotPtr, pluginPath,functionName);
}
void LagrangianDS::setComputeJacobianNNLqFunction(FPtr5 fct)
{
  _pluginJacqNNL->setComputeFunction((void *)fct);
}//computeJacobianNNLqPtr=fct;}
void LagrangianDS::setComputeJacobianNNLqDotFunction(FPtr5 fct)
{
  _pluginJacqDotNNL->setComputeFunction((void *)fct);
}//computeJacobianNNLqDotPtr=fct;}
