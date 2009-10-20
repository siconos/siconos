/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "LagrangianDS.h"
#include "LagrangianDSXML.h"
#include "BlockVector.h"
#include "BlockMatrix.h"

using namespace std;

// Private function to set linked with members of Dynamical top class
void LagrangianDS::connectToDS()
{
  // dim
  _n = 2 * _ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianXF.

  // Initial conditions
  _x0.reset(new BlockVector(_q0, _velocity0));

  // Current State: \f$ x \f$ and rhs = \f$ \dot x \f$

  _x[0].reset(new BlockVector(_q[0], _q[1]));
  _x[1].reset(new BlockVector(_q[1], _q[2]));
  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
}


LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0):
  DynamicalSystem(DS::LNLDS, 2 * newQ0->size()), _ndof(newQ0->size())
{
  zeroPlungin();
  // -- Memory allocation for vector and matrix members --
  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SimpleVector(*_q0));
  _q[1].reset(new SimpleVector(*_velocity0));
  _q[2].reset(new SimpleVector(_ndof));
  mResiduFree.reset(new SimpleVector(getDim()));
  //   mXp.reset(new SimpleVector(getDim()));
  //   mXq.reset(new SimpleVector(getDim()));
  //   mXfree.reset(new SimpleVector(getDim()));
  //   r.reset(new SimpleVector(getDim()));

  // set allocation flags: true for required input, false for others
  _p.resize(3);
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(DS::LNLDS), _ndof(0)
{
  zeroPlungin();
  // Protected constructor - Only call from derived class(es).
  _q.resize(3);
  _p.resize(3);
  // !!! No plug-in connection !!!
}

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(SP::DynamicalSystemXML dsxml):
  DynamicalSystem(dsxml), _ndof(0)
{
  zeroPlungin();
  /*  // -- Lagrangian  xml object --
  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);

  // === Initial conditions ===
  // Warning: ndof is given by q0.size() !!
  if ( ! lgptr->hasQ0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, q0 is a required input");

  q0.reset(new SimpleVector(lgptr->getQ0())); // required node in xml file
  ndof = q0->size();

  if ( ! lgptr->hasVelocity0())
    RuntimeException::selfThrow("LagrangianDS:: xml constructor, v0 is a required input");

  velocity0.reset(new SimpleVector(lgptr->getVelocity0())); // required node in xml file

  if(velocity0->size()!=ndof)
    RuntimeException::selfThrow("LagrangianDS::xml constructor - size of input velocity0 differs from ndof");

  // --- Current State (optional input) ---

  q.resize(3);

  if ( lgptr->hasQ() )
    q[0].reset(new SimpleVector(lgptr->getQ())); // Means that first q is different from q0 ?? Strange case ...
  else
    q[0].reset(new SimpleVector(*q0));           // initialize q with q0
  if ( lgptr->hasVelocity() )
    q[1].reset(new SimpleVector(lgptr->getVelocity0())); // same remark as for q
  else
    q[1].reset(new SimpleVector(*velocity0));

  q[2].reset(new SimpleVector(ndof));
  mResiduFree.reset(new SimpleVector(getDim()));
  p.resize(3);
  // Memories
  if ( lgptr->hasQMemory() ) // qMemory
    qMemory.reset(new SiconosMemory( lgptr->getQMemoryXML() ));

  if ( lgptr->hasVelocityMemory() ) // velocityMemory
    velocityMemory.reset(new SiconosMemory( lgptr->getVelocityMemoryXML() ));

  string plugin;
  // mass
  if( lgptr->isMassPlugin()) // if mass is plugged
    {
      plugin = lgptr->getMassPlugin();
      setComputeMassFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
    }
  else
    mass.reset(new PMMass(lgptr->getMassMatrix()));

  // === Optional inputs ===

  jacobianFInt.resize(2);
  jacobianNNL.resize(2);
  // fInt
  if( lgptr->hasFInt() ) // if fInt is given
    {
      if( lgptr->isFIntPlugin()) // if fInt is plugged
  {
    plugin = lgptr->getFIntPlugin();
    setComputeFIntFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else
  fInt.reset(new PVFint(lgptr->getFIntVector()));
    }

  // _fExt
  if( lgptr->hasFExt() ) // if _fExt is given
    {
      if( lgptr->isFExtPlugin())// if _fExt is plugged
  {
    plugin = lgptr->getFExtPlugin();
    setComputeFExtFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else
  _fExt.reset(new Plugged_Vector_FTime(lgptr->getFExtVector()));
    }

  // NNL
  if( lgptr ->hasNNL())// if NNL is given
    {
      if( lgptr->isNNLPlugin())// if NNL is plugged
  {
    plugin = lgptr->getNNLPlugin();
    setComputeNNLFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else
  NNL.reset(new PVNNL(lgptr->getNNLVector()));
    }

  for(unsigned int i=0;i<2;++i)
    {
      // jacobian(s) of fInt
      if( lgptr ->hasJacobianFInt(i))// if given
  {
    if( lgptr->isJacobianFIntPlugin(i))// if is plugged
      {
        plugin = lgptr->getJacobianFIntPlugin(i);
        setComputeJacobianFIntFunction(i,SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
      }
    else
      jacobianFInt[i].reset(new PMFint(lgptr->getJacobianFIntMatrix(i)));
  }

      // jacobian of NNL
      if( lgptr -> hasJacobianNNL(i)) // if given
  {
    if( lgptr->isJacobianNNLPlugin(i))// if is plugged
      {
        plugin = lgptr->getJacobianNNLPlugin(i);
        setComputeJacobianNNLFunction(i,SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
      }
    else
      jacobianNNL[i].reset(new PMNNL(lgptr->getJacobianNNLMatrix(i)));
  }
  }*/
}

// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum Lagrangian System form: \f$ M\ddot q = p \f$
/**/
LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0, SP::SiconosMatrix newMass):
  DynamicalSystem(DS::LNLDS, 2 * newQ0->size()), _ndof(newQ0->size())
{
  zeroPlungin();
  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  // Mass matrix
  _mass = newMass;

  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SimpleVector(*_q0));
  _q[1].reset(new SimpleVector(*_velocity0));
  _q[2].reset(new SimpleVector(_ndof));
  mResiduFree.reset(new SimpleVector(getDim()));


  _p.resize(3);
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0, const string& massName):
  DynamicalSystem(DS::LNLDS), _ndof(newQ0->size())
{
  zeroPlungin();
  // Initial conditions
  _q0 = newQ0;
  _velocity0 = newVelocity0;

  // Current state
  _q.resize(3);
  _q[0].reset(new SimpleVector(*_q0));
  _q[1].reset(new SimpleVector(*_velocity0));
  _q[2].reset(new SimpleVector(_ndof));
  mResiduFree.reset(new SimpleVector(getDim()));

  // Mass
  setComputeMassFunction(SSL::getPluginName(massName), SSL::getPluginFunctionName(massName));

  _p.resize(3);
}
void LagrangianDS::zeroPlungin()
{
  computeJacobianQNNLPtr = NULL;
  computeJacobianQDotNNLPtr = NULL;
  computeJacobianQFIntPtr = NULL;
  computeJacobianQDotFIntPtr = NULL;
  computeNNLPtr = NULL;
  computeFExtPtr = NULL;
  computeFIntPtr = NULL;
  computeMassPtr = NULL;
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
  if ((_fInt && computeFIntPtr) && (! _jacobianQFInt || ! _jacobianQDotFInt))
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
    output = false;
  }

  // NNL
  if ((_NNL  && computeNNLPtr) && (! _jacobianQNNL || ! _jacobianQDotNNL))
    // ie if NNL is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined NNL but not its Jacobian (according to q and velocity).");
    output = false;
  }

  if (!output) cout << "LagrangianDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

// TEMPORARY FUNCTION: Must be called before this->initialize
void LagrangianDS::initP(const string& simulationType)
{
  if (simulationType == "TimeStepping")
  {
    _p[1].reset(new SimpleVector(_ndof));
    _p[2] = _p[1];
    _p[0] = _p[1];
  }
  else if (simulationType == "EventDriven")
  {
    _p[1].reset(new SimpleVector(_ndof));
    _p[2].reset(new SimpleVector(_ndof));
  }
}

void LagrangianDS::initFL()
{

  _fL.reset(new SimpleVector(_ndof));

  _jacobianQFL.reset(new SimpleMatrix(_ndof, _ndof));
  _jacobianQDotFL.reset(new SimpleMatrix(_ndof, _ndof));
}

void LagrangianDS::initRhs(double time)
{
  _workMatrix.resize(sizeWorkMat);

  // Solve Mq[2]=fL+p.
  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_fL)
  {
    computeFL(time);
    *_q[2] += *_fL;
  }
  computeMass();
  // Copy of Mass into _workMatrix for LU-factorization.

  _workMatrix[invMass].reset(new SimpleMatrix(*_mass));
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);

  bool flag1 = false, flag2 = false;
  if (_jacobianQFL)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianQFL(time);

    _workMatrix[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianQFL));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc10]);
    flag1 = true;
  }

  if (_jacobianQDotFL)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianQDotFL(time);
    _workMatrix[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianQDotFL));
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc11]);
    flag2 = true;
  }

  _workMatrix[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, ZERO));
  _workMatrix[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, IDENTITY));

  if (flag1 && flag2)
    _jacXRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[jacobianXBloc11]));
  else if (flag1) // flag2 = false
    _jacXRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[zeroMatrix]));
  else if (flag2) // flag1 = false
    _jacXRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[jacobianXBloc11]));
  else
    _jacXRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[zeroMatrix]));
}

void LagrangianDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *_q[0] = *_q0;
  *_q[1] = *_velocity0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SimpleVector(1));
  if (computeNNLPtr && !_NNL)
    _NNL.reset(new SimpleVector(_ndof));
  if (computeJacobianQDotNNLPtr && !_jacobianQDotNNL)
    _jacobianQDotNNL.reset(new SimpleMatrix(_ndof, _ndof));
  if (computeJacobianQNNLPtr && ! _jacobianQNNL)
    _jacobianQNNL.reset(new SimpleMatrix(_ndof, _ndof));

  if (computeFExtPtr && !_fExt)
    _fExt.reset(new SimpleVector(_ndof));

  if (computeFIntPtr && ! _fInt)
    _fInt.reset(new SimpleVector(_ndof));
  if (computeJacobianQFIntPtr && !_jacobianQFInt)
    _jacobianQFInt.reset(new SimpleMatrix(_ndof, _ndof));
  if (computeJacobianQDotFIntPtr && !_jacobianQDotFInt)
    _jacobianQDotFInt.reset(new SimpleMatrix(_ndof, _ndof));


  //
  if (!_workFree)
    _workFree.reset(new SimpleVector(getDim()));
  // Memory allocation for fL and its jacobians.
  initFL();

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  if (!_mass)
    _mass.reset(new SimpleMatrix(_ndof, _ndof));
  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  initRhs(time);

}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if (newValue.size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if (! _q[0])
    _q[0].reset(new SimpleVector(newValue));
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
    _q0.reset(new SimpleVector(newValue));
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
    _q[1].reset(new SimpleVector(newValue));
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
    _p[level].reset(new SimpleVector(newValue));
  else
    *(_p[level]) = newValue;
}

void LagrangianDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
{

  if (newPtr->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");
  _p[level] = newPtr;
}
/*
void LagrangianDS::setMass(const PMMass& newValue)
{
  assert(newValue.size(0)==_ndof&&"LagrangianDS - setMass: inconsistent dimensions with problem size for matrix mass.");
  assert(newValue.size(1)==_ndof&&"LagrangianDS - setMass: inconsistent dimensions with problem size for matrix mass.");

  if( !mass  )
    mass.reset(new PMMass(newValue));
  else
    *mass = newValue;
}
void LagrangianDS::setFInt(const PVFint& newValue)
{
  assert(newValue.size()==_ndof&&"LagrangianDS - setFInt: inconsistent dimensions with problem size for input vector fInt");

  if( ! fInt )
    fInt.reset(new PVFint(newValue));
  else
    *_fInt = newValue;
}

void LagrangianDS::setFExt(const SimpleVector& newValue)
{
  assert(newValue.size()==_ndof&&"LagrangianDS - setFExt: inconsistent dimensions with problem size for input vector fExt");

  if( !fExt )
    fExt.reset(new Plugged_Vector_FTime(newValue));
  else
    *_fExt = newValue;
}

void LagrangianDS::setNNL(const PVNNL& newValue)
{
  assert(newValue.size()==_ndof&&"LagrangianDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if( !NNL )
    NNL.reset(new PVNNL(newValue));
  else
    *NNL = newValue;
}

void LagrangianDS::setJacobianFInt(unsigned int i, const PMFint& newValue)
{
  assert(newValue.size(0)==_ndof&&"LagrangianDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");
  assert(newValue.size(1)==_ndof&&"LagrangianDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");

  if( ! jacobianFInt[i] )
    jacobianFInt[i].reset(new PMFint(newValue));
  else
    *_jacobianFInt[i] = newValue;
}

void LagrangianDS::setJacobianNNL(unsigned int i, const PMNNL& newValue)
{
  assert(newValue.size(0)==_ndof&&"LagrangianDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");
  assert(newValue.size(1)==_ndof&&"LagrangianDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");

  if( ! jacobianNNL[i] )
    jacobianNNL[i].reset(new PMNNL(newValue));
  else
    *_jacobianNNL[i] = newValue;
}
*/

void LagrangianDS::computeMass()
{
  if (computeMassPtr)
    (computeMassPtr)(_ndof, &(*_q[0])(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeMass(SP::SiconosVector q2)
{
  if (computeMassPtr)
    (computeMassPtr)(_ndof, &(*q2)(0), &(*_mass)(0, 0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFInt(double time)
{
  if (computeFIntPtr)
    (computeFIntPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeFIntPtr)
    (computeFIntPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeFExt(double time)
{
  if (computeFExtPtr)
    (computeFExtPtr)(time, _ndof, &(*_fExt)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeNNL()
{
  if (computeNNLPtr)
    (computeNNLPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeNNLPtr)
    (computeNNLPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
}

void LagrangianDS::computeJacobianQFInt(double time)
{
  if (computeJacobianQFIntPtr)
    (computeJacobianQFIntPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianQFInt)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianQDotFInt(double time)
{
  if (computeJacobianQDotFIntPtr)
    (computeJacobianQDotFIntPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianQDotFInt)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZFInt(double time){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianQFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQFIntPtr)
    (computeJacobianQFIntPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianQFInt)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianQDotFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQDotFIntPtr)
    (computeJacobianQDotFIntPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianQDotFInt)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianQNNL()
{
  if (computeJacobianQNNLPtr)
    (computeJacobianQNNLPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianQNNL)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianQDotNNL()
{
  if (computeJacobianQDotNNLPtr)
    (computeJacobianQDotNNLPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianQDotNNL)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZNNL(){
//   if(computeJacobianZNNLPtr)
//       (computeJacobianZNNLPtr)(_ndof, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeJacobianQNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQNNLPtr)
    (computeJacobianQNNLPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianQNNL)(0, 0), _z->size(), &(*_z)(0));
}
void LagrangianDS::computeJacobianQDotNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQDotNNLPtr)
    (computeJacobianQDotNNLPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianQDotNNL)(0, 0), _z->size(), &(*_z)(0));
}
// void LagrangianDS::computeJacobianZNNL(unsigned int i, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZNNLPtr)
//     (computeJacobianZNNLPtr)(_ndof, &(*q2)(0), &(*velocity2)(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void LagrangianDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  *_q[2] = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_fL)
  {
    computeFL(time);
    *_q[2] += *_fL;
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

void LagrangianDS::computeJacobianXRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  if (!isDSup)
    computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children
  *_workMatrix[invMass] = *_mass;

  if (_jacobianQFL)
  {
    SP::SiconosMatrix bloc10 = _jacXRhs->block(1, 0);
    computeJacobianQFL(time);
    *bloc10 = *_jacobianQFL;
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc10);
  }

  if (_jacobianQDotFL)
  {
    SP::SiconosMatrix bloc11 = _jacXRhs->block(1, 1);
    computeJacobianQDotFL(time);
    *bloc11 = *_jacobianQDotFL;
    _workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeFL(double time)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_fL)
  {
    // 1 - Computes the required functions
    computeFInt(time);
    computeFExt(time);
    computeNNL();

    // 2 - set fL = fExt - fInt - NNL

    // seems ok.
    if (_fL.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      _fL->zero();

      if (_fInt)
        *_fL -= *_fInt;

      if (_fExt)
        *_fL += *_fExt;

      if (_NNL)
        *_fL -= *_NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeFL(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_fL)
  {
    // 1 - Computes the required functions
    computeFInt(time, q2, v2);
    computeFExt(time);
    computeNNL(q2, v2);

    // seems ok.
    if (_fL.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      _fL->zero();

      if (_fInt)
        *_fL -= *_fInt;

      if (_fExt)
        *_fL += *_fExt;

      if (_NNL)
        *_fL -= *_NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeJacobianQFL(double time)
{
  if (_jacobianQFL)
  {
    computeJacobianQFInt(time);
    computeJacobianQNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianQFL->zero();
      if (_jacobianQFInt)
        *_jacobianQFL -= *_jacobianQFInt;
      if (_jacobianQNNL)
        *_jacobianQFL -= *_jacobianQNNL;
    }
  }
  //else nothing.
}
void LagrangianDS::computeJacobianQDotFL(double time)
{
  if (_jacobianQDotFL)
  {
    computeJacobianQDotFInt(time);
    computeJacobianQDotNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianQDotFL->zero();
      if (_jacobianQDotFInt)
        *_jacobianQDotFL -= *_jacobianQDotFInt;
      if (_jacobianQDotNNL)
        *_jacobianQDotFL -= *_jacobianQDotNNL;
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
  /*  if(!dsxml)
    RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DynamicalSystemXML does not exist");

  SP::LagrangianDSXML lgptr = boost::static_pointer_cast <LagrangianDSXML>(dsxml);
  lgptr->setMassPlugin(mass->getPluginName());
  lgptr->setQ( *_q[0] );
  lgptr->setQ0( *_q0 );
  lgptr->setQMemory( *_qMemory );
  lgptr->setVelocity( *q[1] );
  lgptr->setVelocity0( *velocity0 );
  lgptr->setVelocityMemory( *velocityMemory );

  // FExt
  if( lgptr->hasFExt() )
    {
      if( !lgptr->isFExtPlugin())
  lgptr->setFExtVector( *_fExt );
    }
  else
    lgptr->setFExtPlugin(fExt->getPluginName());

  // FInt
  if( lgptr->hasFInt() )
    {
      if( !lgptr->isFIntPlugin())
  if( fInt->size() > 0 )
    lgptr->setFIntVector(*_fInt);
  else cout<<"Warning : FInt can't be saved, the FInt vector is not defined."<<endl;
    }
  else
    lgptr->setFIntPlugin(fInt->getPluginName());

  for(unsigned int i=0;i<2;++i)
    {
      // Jacobian FInt
      if( lgptr->hasJacobianFInt(i) )
  {
    if( !lgptr->isJacobianFIntPlugin(i))
      lgptr->setJacobianFIntMatrix(i,*_jacobianFInt[i]);
  }
      else
  lgptr->setJacobianFIntPlugin(i,jacobianFInt[i]->getPluginName());

      // JacobianQNNL
      if( lgptr->hasJacobianNNL(i) )
  {
    if( !lgptr->isJacobianNNLPlugin(i))
      lgptr->setJacobianNNLMatrix(i, *_jacobianNNL[i] );
  }
      else
  lgptr->setJacobianNNLPlugin(i,jacobianNNL[i]->getPluginName());
    }
  // NNL
  if( lgptr->hasNNL() )
    {
      if( !lgptr->isNNLPlugin())
  lgptr->setNNLVector(*_NNL);
    }
  else
  lgptr->setNNLPlugin(NNL->getPluginName());*/
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
  cout << "- p " << endl;
  if (_p[2]) _p[2]->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {
    _qMemory.reset(new SiconosMemory(steps));
    _velocityMemory.reset(new SiconosMemory(steps));
    swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{

  _xMemory->swap(_x[0]);
  _qMemory->swap(_q[0]);
  _velocityMemory->swap(_q[1]);
  // initialization of the reaction force due to the non smooth law
  _p[1]->zero();
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

void LagrangianDS::computeQFree(double time, unsigned int level, SP::SiconosVector qFreeOut)
{
  // to compute qFree, derivative number level. Mainly used in EventDriven to compute second derivative
  // of q for Output y computation.

  if (qFreeOut->size() != _ndof)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Wrong size for output (different from _ndof)");


  if (level != 2)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Only implemented for second derivative.");

  // Warning: we suppose that all other operators are up to date (FInt, FExt ...)

  qFreeOut->zero();
  if (_fInt)
    *qFreeOut -= *_fInt;
  if (_fExt)
    *qFreeOut += *_fExt;
  if (_NNL)
    *qFreeOut -= *_NNL;

  _workMatrix[invMass]->PLUForwardBackwardInPlace(*qFreeOut);
}

void LagrangianDS::resetNonSmoothPart()
{
  _p[1]->zero();
}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_p[1]);
  *_q[1] += *_p[1];  // v+ = v- + p
}

void LagrangianDS::setComputeNNLFunction(const std::string& pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeNNLPtr, pluginPath, functionName);
}
void LagrangianDS::setComputeNNLFunction(FPtr5 fct)
{
  computeNNLPtr = fct;
}
void LagrangianDS::setComputeJacobianQFIntFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianQDotFIntFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianQDotFIntPtr, pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianQFIntFunction(FPtr6 fct)
{
  computeJacobianQFIntPtr = fct;
}
void LagrangianDS::setComputeJacobianQDotFIntFunction(FPtr6 fct)
{
  computeJacobianQDotFIntPtr = fct;
}
void LagrangianDS::setComputeJacobianQNNLFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianQNNLPtr, pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianQDotNNLFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianQDotNNLPtr, pluginPath, functionName);
}
void LagrangianDS::setComputeJacobianQNNLFunction(FPtr5 fct)
{
  computeJacobianQNNLPtr = fct;
}
void LagrangianDS::setComputeJacobianQDotNNLFunction(FPtr5 fct)
{
  computeJacobianQDotNNLPtr = fct;
}
