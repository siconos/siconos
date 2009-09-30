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
  n = 2 * ndof;

  // All links between DS and LagrangianDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianXF.

  // Initial conditions
  x0.reset(new BlockVector(q0, velocity0));

  // Current State: \f$ x \f$ and rhs = \f$ \dot x \f$

  x[0].reset(new BlockVector(q[0], q[1]));
  x[1].reset(new BlockVector(q[1], q[2]));
  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
}


LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0):
  DynamicalSystem(DS::LNLDS, 2 * newQ0->size()), ndof(newQ0->size())
{
  zeroPlungin();
  // -- Memory allocation for vector and matrix members --
  // Initial conditions
  q0 = newQ0;
  velocity0 = newVelocity0;

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));
  mResiduFree.reset(new SimpleVector(getDim()));
  //   mXp.reset(new SimpleVector(getDim()));
  //   mXq.reset(new SimpleVector(getDim()));
  //   mXfree.reset(new SimpleVector(getDim()));
  //   r.reset(new SimpleVector(getDim()));

  // set allocation flags: true for required input, false for others
  p.resize(3);
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(DS::LNLDS), ndof(0)
{
  zeroPlungin();
  // Protected constructor - Only call from derived class(es).
  q.resize(3);
  p.resize(3);
  // !!! No plug-in connection !!!
}

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(SP::DynamicalSystemXML dsxml):
  DynamicalSystem(dsxml), ndof(0)
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

  // fExt
  if( lgptr->hasFExt() ) // if fExt is given
    {
      if( lgptr->isFExtPlugin())// if fExt is plugged
  {
    plugin = lgptr->getFExtPlugin();
    setComputeFExtFunction(SSL::getPluginName( plugin ), SSL::getPluginFunctionName( plugin ));
  }
      else
  fExt.reset(new Plugged_Vector_FTime(lgptr->getFExtVector()));
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
  DynamicalSystem(DS::LNLDS, 2 * newQ0->size()), ndof(newQ0->size())
{
  zeroPlungin();
  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  // Mass matrix
  mass = newMass;

  // Initial conditions
  q0 = newQ0;
  velocity0 = newVelocity0;

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));
  mResiduFree.reset(new SimpleVector(getDim()));


  p.resize(3);
}

// From a set of data - Mass loaded from a plugin
// This constructor leads to the minimum Lagrangian System form: \f$ M(q)\ddot q = p \f$
LagrangianDS::LagrangianDS(SP::SiconosVector newQ0, SP::SiconosVector newVelocity0, const string& massName):
  DynamicalSystem(DS::LNLDS), ndof(newQ0->size())
{
  zeroPlungin();
  // Initial conditions
  q0 = newQ0;
  velocity0 = newVelocity0;

  // Current state
  q.resize(3);
  q[0].reset(new SimpleVector(*q0));
  q[1].reset(new SimpleVector(*velocity0));
  q[2].reset(new SimpleVector(ndof));
  mResiduFree.reset(new SimpleVector(getDim()));

  // Mass
  setComputeMassFunction(SSL::getPluginName(massName), SSL::getPluginFunctionName(massName));

  p.resize(3);
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
  if (ndof == 0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0
  if (! q0 || ! velocity0)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }

  // Mass
  if (! mass)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - Mass not set.");
    output = false;
  }

  // fInt
  if ((fInt && computeFIntPtr) && (! jacobianQFInt || ! jacobianQDotFInt))
    // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("LagrangianDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
    output = false;
  }

  // NNL
  if ((NNL  && computeNNLPtr) && (! jacobianQNNL || ! jacobianQDotNNL))
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
    p[1].reset(new SimpleVector(ndof));
    p[2] = p[1];
    p[0] = p[1];
  }
  else if (simulationType == "EventDriven")
  {
    p[1].reset(new SimpleVector(ndof));
    p[2].reset(new SimpleVector(ndof));
  }
}

void LagrangianDS::initFL()
{

  fL.reset(new SimpleVector(ndof));

  jacobianQFL.reset(new SimpleMatrix(ndof, ndof));
  jacobianQDotFL.reset(new SimpleMatrix(ndof, ndof));
}

void LagrangianDS::initRhs(double time)
{
  workMatrix.resize(sizeWorkMat);

  // Solve Mq[2]=fL+p.
  *q[2] = *(p[2]); // Warning: r/p update is done in Interactions/Relations

  if (fL)
  {
    computeFL(time);
    *q[2] += *fL;
  }
  computeMass();
  // Copy of Mass into workMatrix for LU-factorization.

  workMatrix[invMass].reset(new SimpleMatrix(*mass));
  workMatrix[invMass]->PLUForwardBackwardInPlace(*q[2]);

  bool flag1 = false, flag2 = false;
  if (jacobianQFL)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianQFL(time);

    workMatrix[jacobianXBloc10].reset(new SimpleMatrix(*jacobianQFL));
    workMatrix[invMass]->PLUForwardBackwardInPlace(*workMatrix[jacobianXBloc10]);
    flag1 = true;
  }

  if (jacobianQDotFL)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianQDotFL(time);
    workMatrix[jacobianXBloc11].reset(new SimpleMatrix(*jacobianQDotFL));
    workMatrix[invMass]->PLUForwardBackwardInPlace(*workMatrix[jacobianXBloc11]);
    flag2 = true;
  }

  workMatrix[zeroMatrix].reset(new SimpleMatrix(ndof, ndof, ZERO));
  workMatrix[idMatrix].reset(new SimpleMatrix(ndof, ndof, IDENTITY));

  if (flag1 && flag2)
    jacobianXRhs.reset(new BlockMatrix(workMatrix[zeroMatrix], workMatrix[idMatrix], workMatrix[jacobianXBloc10], workMatrix[jacobianXBloc11]));
  else if (flag1) // flag2 = false
    jacobianXRhs.reset(new BlockMatrix(workMatrix[zeroMatrix], workMatrix[idMatrix], workMatrix[jacobianXBloc10], workMatrix[zeroMatrix]));
  else if (flag2) // flag1 = false
    jacobianXRhs.reset(new BlockMatrix(workMatrix[zeroMatrix], workMatrix[idMatrix], workMatrix[zeroMatrix], workMatrix[jacobianXBloc11]));
  else
    jacobianXRhs.reset(new BlockMatrix(workMatrix[zeroMatrix], workMatrix[idMatrix], workMatrix[zeroMatrix], workMatrix[zeroMatrix]));
}

void LagrangianDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *q[0] = *q0;
  *q[1] = *velocity0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! z)
    z.reset(new SimpleVector(1));
  if (computeNNLPtr && !NNL)
    NNL.reset(new SimpleVector(ndof));
  if (computeJacobianQDotNNLPtr && !jacobianQDotNNL)
    jacobianQDotNNL.reset(new SimpleMatrix(ndof, ndof));
  if (computeJacobianQNNLPtr && ! jacobianQNNL)
    jacobianQNNL.reset(new SimpleMatrix(ndof, ndof));

  if (computeFExtPtr && !fExt)
    fExt.reset(new SimpleVector(ndof));

  if (computeFIntPtr && ! fInt)
    fInt.reset(new SimpleVector(ndof));
  if (computeJacobianQFIntPtr && !jacobianQFInt)
    jacobianQFInt.reset(new SimpleMatrix(ndof, ndof));
  if (computeJacobianQDotFIntPtr && !jacobianQDotFInt)
    jacobianQDotFInt.reset(new SimpleMatrix(ndof, ndof));



  // Memory allocation for fL and its jacobians.
  initFL();

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  if (!mass)
    mass.reset(new SimpleMatrix(ndof, ndof));
  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  initRhs(time);

}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQ(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ: inconsistent input vector size ");

  if (! q[0])
    q[0].reset(new SimpleVector(newValue));
  else
    *q[0] = newValue;
}

void LagrangianDS::setQPtr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQPtr: inconsistent input vector size ");
  q[0] = newPtr;

}

void LagrangianDS::setQ0(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0: inconsistent input vector size ");

  if (! q0)
    q0.reset(new SimpleVector(newValue));
  else
    *q0 = newValue;
}

void LagrangianDS::setQ0Ptr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setQ0Ptr: inconsistent input vector size ");
  q0 = newPtr;
}

void LagrangianDS::setVelocity(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity: inconsistent input vector size ");

  if (! q[1])
    q[1].reset(new SimpleVector(newValue));
  else
    *q[1] = newValue;
}

void LagrangianDS::setVelocityPtr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocityPtr: inconsistent input vector size ");
  q[1] = newPtr;
}

void LagrangianDS::setVelocity0(const SiconosVector& newValue)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0: inconsistent input vector size ");

  if (! velocity0)
    velocity0.reset(new SimpleVector(newValue));
  else
    *velocity0 = newValue;
}

void LagrangianDS::setVelocity0Ptr(SP::SiconosVector newPtr)
{
  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setVelocity0Ptr: inconsistent input vector size ");
  velocity0 = newPtr;
}

SP::SiconosVector LagrangianDS::getAccelerationPtr() const
{
  return q[2];
}

void LagrangianDS::setP(const SiconosVector& newValue, unsigned int level)
{
  if (newValue.size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setP: inconsistent input vector size ");

  if (! p[level])
    p[level].reset(new SimpleVector(newValue));
  else
    *(p[level]) = newValue;
}

void LagrangianDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
{

  if (newPtr->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS - setPPtr: inconsistent input vector size ");
  p[level] = newPtr;
}
/*
void LagrangianDS::setMass(const PMMass& newValue)
{
  assert(newValue.size(0)==ndof&&"LagrangianDS - setMass: inconsistent dimensions with problem size for matrix mass.");
  assert(newValue.size(1)==ndof&&"LagrangianDS - setMass: inconsistent dimensions with problem size for matrix mass.");

  if( !mass  )
    mass.reset(new PMMass(newValue));
  else
    *mass = newValue;
}
void LagrangianDS::setFInt(const PVFint& newValue)
{
  assert(newValue.size()==ndof&&"LagrangianDS - setFInt: inconsistent dimensions with problem size for input vector fInt");

  if( ! fInt )
    fInt.reset(new PVFint(newValue));
  else
    *fInt = newValue;
}

void LagrangianDS::setFExt(const SimpleVector& newValue)
{
  assert(newValue.size()==ndof&&"LagrangianDS - setFExt: inconsistent dimensions with problem size for input vector fExt");

  if( !fExt )
    fExt.reset(new Plugged_Vector_FTime(newValue));
  else
    *fExt = newValue;
}

void LagrangianDS::setNNL(const PVNNL& newValue)
{
  assert(newValue.size()==ndof&&"LagrangianDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if( !NNL )
    NNL.reset(new PVNNL(newValue));
  else
    *NNL = newValue;
}

void LagrangianDS::setJacobianFInt(unsigned int i, const PMFint& newValue)
{
  assert(newValue.size(0)==ndof&&"LagrangianDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");
  assert(newValue.size(1)==ndof&&"LagrangianDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");

  if( ! jacobianFInt[i] )
    jacobianFInt[i].reset(new PMFint(newValue));
  else
    *jacobianFInt[i] = newValue;
}

void LagrangianDS::setJacobianNNL(unsigned int i, const PMNNL& newValue)
{
  assert(newValue.size(0)==ndof&&"LagrangianDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");
  assert(newValue.size(1)==ndof&&"LagrangianDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");

  if( ! jacobianNNL[i] )
    jacobianNNL[i].reset(new PMNNL(newValue));
  else
    *jacobianNNL[i] = newValue;
}
*/

void LagrangianDS::computeMass()
{
  if (computeMassPtr)
    (computeMassPtr)(ndof, &(*q[0])(0), &(*mass)(0, 0), z->size(), &(*z)(0));
}

void LagrangianDS::computeMass(SP::SiconosVector q2)
{
  if (computeMassPtr)
    (computeMassPtr)(ndof, &(*q2)(0), &(*mass)(0, 0), z->size(), &(*z)(0));
}

void LagrangianDS::computeFInt(double time)
{
  if (computeFIntPtr)
    (computeFIntPtr)(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*fInt)(0), z->size(), &(*z)(0));
}
void LagrangianDS::computeFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeFIntPtr)
    (computeFIntPtr)(time, ndof, &(*q2)(0), &(*velocity2)(0), &(*fInt)(0), z->size(), &(*z)(0));
}

void LagrangianDS::computeFExt(double time)
{
  if (computeFExtPtr)
    (computeFExtPtr)(time, ndof, &(*fExt)(0), z->size(), &(*z)(0));
}

void LagrangianDS::computeNNL()
{
  if (computeNNLPtr)
    (computeNNLPtr)(ndof, &(*q[0])(0), &(*q[1])(0), &(*NNL)(0), z->size(), &(*z)(0));
}

void LagrangianDS::computeNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeNNLPtr)
    (computeNNLPtr)(ndof, &(*q2)(0), &(*velocity2)(0), &(*NNL)(0), z->size(), &(*z)(0));
}

void LagrangianDS::computeJacobianQFInt(double time)
{
  if (computeJacobianQFIntPtr)
    (computeJacobianQFIntPtr)(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianQFInt)(0, 0), z->size(), &(*z)(0));
}
void LagrangianDS::computeJacobianQDotFInt(double time)
{
  if (computeJacobianQDotFIntPtr)
    (computeJacobianQDotFIntPtr)(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianQDotFInt)(0, 0), z->size(), &(*z)(0));
}
// void LagrangianDS::computeJacobianZFInt(double time){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianFInt[i])(0,0), z->size(), &(*z)(0));
// }

void LagrangianDS::computeJacobianQFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQFIntPtr)
    (computeJacobianQFIntPtr)(time, ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQFInt)(0, 0), z->size(), &(*z)(0));
}
void LagrangianDS::computeJacobianQDotFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQDotFIntPtr)
    (computeJacobianQDotFIntPtr)(time, ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQDotFInt)(0, 0), z->size(), &(*z)(0));
}
// void LagrangianDS::computeJacobianZFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianFInt[i])(0,0), z->size(), &(*z)(0));
// }

void LagrangianDS::computeJacobianQNNL()
{
  if (computeJacobianQNNLPtr)
    (computeJacobianQNNLPtr)(ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianQNNL)(0, 0), z->size(), &(*z)(0));
}
void LagrangianDS::computeJacobianQDotNNL()
{
  if (computeJacobianQDotNNLPtr)
    (computeJacobianQDotNNLPtr)(ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianQDotNNL)(0, 0), z->size(), &(*z)(0));
}
// void LagrangianDS::computeJacobianZNNL(){
//   if(computeJacobianZNNLPtr)
//       (computeJacobianZNNLPtr)(ndof, &(*q[0])(0), &(*q[1])(0), &(*jacobianNNL[i])(0,0), z->size(), &(*z)(0));
// }

void LagrangianDS::computeJacobianQNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQNNLPtr)
    (computeJacobianQNNLPtr)(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQNNL)(0, 0), z->size(), &(*z)(0));
}
void LagrangianDS::computeJacobianQDotNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianQDotNNLPtr)
    (computeJacobianQDotNNLPtr)(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianQDotNNL)(0, 0), z->size(), &(*z)(0));
}
// void LagrangianDS::computeJacobianZNNL(unsigned int i, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZNNLPtr)
//     (computeJacobianZNNLPtr)(ndof, &(*q2)(0), &(*velocity2)(0), &(*jacobianNNL[i])(0,0), z->size(), &(*z)(0));
// }

void LagrangianDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  *q[2] = *(p[2]); // Warning: r/p update is done in Interactions/Relations

  if (fL)
  {
    computeFL(time);
    *q[2] += *fL;
  }

  // mass and inv(mass) computatiton
  if (!isDSup) // if it is necessary to re-compute mass, FInt ..., ie if they have not been compiled during the present time step
    computeMass();


  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into workMatrix[invMass].
  // -- Case 2: mass is not constant, we copy it into workMatrix[invMass]
  // Then we proceed with PLUForwardBackward.

  //  if(mass->isPlugged()) : mass may be not plugged in LagrangianDS children
  *workMatrix[invMass] = *mass;

  workMatrix[invMass]->PLUForwardBackwardInPlace(*q[2]);

}

void LagrangianDS::computeJacobianXRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  if (!isDSup)
    computeMass();

  //  if(mass->isPlugged()) : mass may b not plugged in LagrangianDS children
  *workMatrix[invMass] = *mass;

  if (jacobianQFL)
  {
    SP::SiconosMatrix bloc10 = jacobianXRhs->getBlockPtr(1, 0);
    computeJacobianQFL(time);
    *bloc10 = *jacobianQFL;
    workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc10);
  }

  if (jacobianQDotFL)
  {
    SP::SiconosMatrix bloc11 = jacobianXRhs->getBlockPtr(1, 1);
    computeJacobianQDotFL(time);
    *bloc11 = *jacobianQDotFL;
    workMatrix[invMass]->PLUForwardBackwardInPlace(*bloc11);
  }
}

void LagrangianDS::computeFL(double time)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (fL)
  {
    // 1 - Computes the required functions
    computeFInt(time);
    computeFExt(time);
    computeNNL();

    // 2 - set fL = fExt - fInt - NNL

    // seems ok.
    if (fL.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      fL->zero();

      if (fInt)
        *fL -= *fInt;

      if (fExt)
        *fL += *fExt;

      if (NNL)
        *fL -= *NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeFL(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (fL)
  {
    // 1 - Computes the required functions
    computeFInt(time, q2, v2);
    computeFExt(time);
    computeNNL(q2, v2);

    // seems ok.
    if (fL.use_count() == 1)
    {
      //if not that means that fL is already (pointer-)connected with
      // either fInt, NNL OR fExt.
      fL->zero();

      if (fInt)
        *fL -= *fInt;

      if (fExt)
        *fL += *fExt;

      if (NNL)
        *fL -= *NNL;
    }
  }
  // else nothing.
}

void LagrangianDS::computeJacobianQFL(double time)
{
  if (jacobianQFL)
  {
    computeJacobianQFInt(time);
    computeJacobianQNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      jacobianQFL->zero();
      if (jacobianQFInt)
        *jacobianQFL -= *jacobianQFInt;
      if (jacobianQNNL)
        *jacobianQFL -= *jacobianQNNL;
    }
  }
  //else nothing.
}
void LagrangianDS::computeJacobianQDotFL(double time)
{
  if (jacobianQDotFL)
  {
    computeJacobianQDotFInt(time);
    computeJacobianQDotNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      jacobianQDotFL->zero();
      if (jacobianQDotFInt)
        *jacobianQDotFL -= *jacobianQDotFInt;
      if (jacobianQDotNNL)
        *jacobianQDotFL -= *jacobianQDotNNL;
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
  lgptr->setQ( *q[0] );
  lgptr->setQ0( *q0 );
  lgptr->setQMemory( *qMemory );
  lgptr->setVelocity( *q[1] );
  lgptr->setVelocity0( *velocity0 );
  lgptr->setVelocityMemory( *velocityMemory );

  // FExt
  if( lgptr->hasFExt() )
    {
      if( !lgptr->isFExtPlugin())
  lgptr->setFExtVector( *fExt );
    }
  else
    lgptr->setFExtPlugin(fExt->getPluginName());

  // FInt
  if( lgptr->hasFInt() )
    {
      if( !lgptr->isFIntPlugin())
  if( fInt->size() > 0 )
    lgptr->setFIntVector(*fInt);
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
      lgptr->setJacobianFIntMatrix(i,*jacobianFInt[i]);
  }
      else
  lgptr->setJacobianFIntPlugin(i,jacobianFInt[i]->getPluginName());

      // JacobianQNNL
      if( lgptr->hasJacobianNNL(i) )
  {
    if( !lgptr->isJacobianNNLPlugin(i))
      lgptr->setJacobianNNLMatrix(i, *jacobianNNL[i] );
  }
      else
  lgptr->setJacobianNNLPlugin(i,jacobianNNL[i]->getPluginName());
    }
  // NNL
  if( lgptr->hasNNL() )
    {
      if( !lgptr->isNNLPlugin())
  lgptr->setNNLVector(*NNL);
    }
  else
  lgptr->setNNLPlugin(NNL->getPluginName());*/
}

void LagrangianDS::display() const
{
  cout << "=====> Lagrangian System display (number: " << number << ")." << endl;
  cout << "- ndof : " << ndof << endl;
  cout << "- q " << endl;
  if (q[0]) q[0]->display();
  else cout << "-> NULL" << endl;
  cout << "- q0 " << endl;
  if (q0) q0->display();
  cout << "- v " << endl;
  if (q[1]) q[1]->display();
  else cout << "-> NULL" << endl;
  cout << "- v0 " << endl;
  if (velocity0) velocity0->display();
  else cout << "-> NULL" << endl;
  cout << "- p " << endl;
  if (p[2]) p[2]->display();
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
    qMemory.reset(new SiconosMemory(steps));
    velocityMemory.reset(new SiconosMemory(steps));
    swapInMemory();
  }
}

void LagrangianDS::swapInMemory()
{

  xMemory->swap(x[0]);
  qMemory->swap(q[0]);
  velocityMemory->swap(q[1]);
  // initialization of the reaction force due to the non smooth law
  p[1]->zero();
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

  if (qFreeOut->size() != ndof)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Wrong size for output (different from ndof)");


  if (level != 2)
    RuntimeException::selfThrow("LagrangianDS::computeQFree - Only implemented for second derivative.");

  // Warning: we suppose that all other operators are up to date (FInt, FExt ...)

  qFreeOut->zero();
  if (fInt)
    *qFreeOut -= *fInt;
  if (fExt)
    *qFreeOut += *fExt;
  if (NNL)
    *qFreeOut -= *NNL;

  workMatrix[invMass]->PLUForwardBackwardInPlace(*qFreeOut);
}

void LagrangianDS::resetNonSmoothPart()
{
  p[1]->zero();
}

void LagrangianDS::computePostImpactVelocity()
{
  // When this function is call, q[1] is supposed to be pre-impact velocity.
  // We solve M(v+ - v-) = p - The result is saved in(place of) p[1].
  workMatrix[invMass]->PLUForwardBackwardInPlace(*p[1]);
  *q[1] += *p[1];  // v+ = v- + p
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
