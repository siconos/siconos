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
#include "NewtonEulerDS.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"

using namespace std;

// Private function to set linked with members of Dynamical top class
void NewtonEulerDS::connectToDS()
{
  // dim
  _n = 2 * _ndof;

  // All links between DS and NewtonEulerDS class members are pointer links, which means
  // that no useless memory is allocated when connection is established.
  // One exception: zero and identity matrices, used to filled in M and jacobianXF.

  // Initial conditions
  //  _x0.reset(new BlockVector(_q0,_v0));

  // Current State: \f$ x \f$ and rhs = \f$ \dot x \f$

  //  _x[0].reset(new BlockVector(_q[0],_q[1]));
  //  _x[1].reset(new BlockVector(_q[1],_q[2]));
  // Everything concerning rhs and its jacobian is handled in initRhs and computeXXX related functions.
}




// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum NewtonEuler System form: \f$ M\ddot q = p \f$
/*
Q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
Velocity0 : contains the initial velocity of center of mass and the omega initial. (dim(Velocity0)=6)
*/
NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix):
  DynamicalSystem(DS::NENLDS, 6)
{
  _p.resize(3);
  zeroPlugin();
  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  _mass = mass;
  _ndof = 3;
  _qDim = 7;
  _n = 6;

  // Initial conditions
  _q0 = Q0;
  _v0 = Velocity0;

  // Current state
  _q.reset(new SimpleVector(_qDim));
  _v.reset(new SimpleVector(_n));
  (*_q) = (*_q0);
  _dotq.reset(new SimpleVector(_qDim));
  _residuFree.reset(new SimpleVector(_n));
  _W.reset(new SimpleMatrix(_n, _n));
  _luW.reset(new SimpleMatrix(_n, _n));
  _W->zero();
  _W->setValue(0, 0, _mass);
  _W->setValue(1, 1, _mass);
  _W->setValue(2, 2, _mass);
  _I = inertialMatrix;
  Index dimIndex(2);
  dimIndex[0] = 3;
  dimIndex[1] = 3;
  Index startIndex(4);
  startIndex[0] = 0;
  startIndex[1] = 0;
  startIndex[2] = 3;
  startIndex[3] = 3;
  setBlock(_I, _W, dimIndex, startIndex);
  *_luW = *_W;


  _T.reset(new SimpleMatrix(_qDim, _n));
  _T->zero();
  _T->setValue(0, 0, 1.0);
  _T->setValue(1, 1, 1.0);
  _T->setValue(2, 2, 1.0);
  updateT();
}
void NewtonEulerDS::zeroPlugin()
{
  computeJacobianqNNLPtr = NULL;
  computeJacobianqDotNNLPtr = NULL;
  computeJacobianqFIntPtr = NULL;
  computeJacobianqDotFIntPtr = NULL;
  computeNNLPtr = NULL;
  computeFExtPtr = NULL;
  //  computeFIntPtr=NULL;
}

// Destructor
NewtonEulerDS::~NewtonEulerDS()
{
}

bool NewtonEulerDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof
  if (_ndof == 0)
  {
    RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - number of degrees of freedom is equal to 0.");
    output = false;
  }

  // q0 and velocity0
  if (! _q0 || ! _v0)
  {
    RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }


  // fInt
  //   if( ( _fInt && computeFIntPtr) && ( ! _jacobianqFInt || ! _jacobianqDotFInt ) )
  //     // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  //     {
  //       RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
  //       output = false;
  //     }

  // NNL
  if ((_NNL  && computeNNLPtr) && (! _jacobianqNNL || ! _jacobianqDotNNL))
    // ie if NNL is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  {
    RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - You defined NNL but not its Jacobian (according to q and velocity).");
    output = false;
  }

  if (!output) cout << "NewtonEulerDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

// TEMPORARY FUNCTION: Must be called before this->initialize
void NewtonEulerDS::initP(const string& simulationType)
{
  if (simulationType == "TimeStepping")
  {
    _p[1].reset(new SimpleVector(_n));
    _p[2] = _p[1];
    _p[0] = _p[1];
  }
  else if (simulationType == "EventDriven")
  {
    _p[1].reset(new SimpleVector(_n));
    _p[2].reset(new SimpleVector(_n));
  }
}

void NewtonEulerDS::initFL()
{

  _fL.reset(new SimpleVector(_n));

  _jacobianqFL.reset(new SimpleMatrix(_n, _qDim));
  _jacobianqDotFL.reset(new SimpleMatrix(_n, _qDim));
}

void NewtonEulerDS::initRhs(double time)
{
  //  _workMatrix.resize(sizeWorkMat);

  // Solve Mq[2]=fL+p.
  //*_q = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_fL)
  {
    computeFL(time);
    //      *_q += *_fL;
  }

  // Copy of Mass into _workMatrix for LU-factorization.

  bool flag1 = false, flag2 = false;
  if (_jacobianqFL)
  {
    // Solve MjacobianX(1,0) = jacobianFL[0]
    computeJacobianqFL(time);

    //      _workMatrix[jacobianXBloc10].reset(new SimpleMatrix(*_jacobianqFL));
    flag1 = true;
  }

  if (_jacobianqDotFL)
  {
    // Solve MjacobianX(1,1) = jacobianFL[1]
    computeJacobianqDotFL(time);
    //      _workMatrix[jacobianXBloc11].reset(new SimpleMatrix(*_jacobianqDotFL));
    //      _workMatrix[invMass]->PLUForwardBackwardInPlace(*_workMatrix[jacobianXBloc11]);
    flag2 = true;
  }

  //   _workMatrix[zeroMatrix].reset(new SimpleMatrix(_ndof, _ndof, ZERO));
  //   _workMatrix[idMatrix].reset(new SimpleMatrix(_ndof, _ndof, IDENTITY));

  //   if(flag1&&flag2)
  //     _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[jacobianXBloc11]));
  //   else if(flag1) // flag2 = false
  //     _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[jacobianXBloc10], _workMatrix[zeroMatrix]));
  //   else if(flag2) // flag1 = false
  //     _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[jacobianXBloc11]));
  //   else
  //     _jacxRhs.reset(new BlockMatrix(_workMatrix[zeroMatrix], _workMatrix[idMatrix], _workMatrix[zeroMatrix], _workMatrix[zeroMatrix]));
}

void NewtonEulerDS::initialize(const string& simulationType, double time, unsigned int sizeOfMemory)
{
  // Memory allocation for p[0], p[1], p[2].
  initP(simulationType);

  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *_q = *_q0;
  *_v = *_v0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SimpleVector(1));
  if (computeNNLPtr && !_NNL)
    _NNL.reset(new SimpleVector(_n));
  if (computeJacobianqDotNNLPtr && !_jacobianqDotNNL)
    _jacobianqDotNNL.reset(new SimpleMatrix(_n, _n));
  if (computeJacobianqNNLPtr && ! _jacobianqNNL)
    _jacobianqNNL.reset(new SimpleMatrix(_n, _n));

  if (computeFExtPtr && !_fExt)
    _fExt.reset(new SimpleVector(_n));

  //   if ( computeFIntPtr && ! _fInt)
  //     _fInt.reset(new SimpleVector(_n));
  //   if (computeJacobianqFIntPtr && !_jacobianqFInt)
  //     _jacobianqFInt.reset(new SimpleMatrix(_n,_n));
  //   if (computeJacobianqDotFIntPtr && !_jacobianqDotFInt)
  //     _jacobianqDotFInt.reset(new SimpleMatrix(_n,_n));


  //
  if (!_workFree)
    _workFree.reset(new SimpleVector(getDim()));
  // Memory allocation for fL and its jacobians.
  initFL();

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  initRhs(time);

}

// --- GETTERS/SETTERS ---

// void NewtonEulerDS::setQ(const SiconosVector& newValue)
// {
//   if(newValue.size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setQ: inconsistent input vector size ");

//   if( ! _q[0] )
//     _q[0].reset(new SimpleVector(newValue));
//   else
//     *_q[0] = newValue;
// }

// void NewtonEulerDS::setQPtr(SP::SiconosVector newPtr)
// {
//   if(newPtr->size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setQPtr: inconsistent input vector size ");
//   _q[0] = newPtr;

// }

// void NewtonEulerDS::setQ0(const SiconosVector& newValue)
// {
//   if(newValue.size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setQ0: inconsistent input vector size ");

//   if( ! _q0 )
//     _q0.reset(new SimpleVector(newValue));
//   else
//     *_q0 = newValue;
// }

// void NewtonEulerDS::setQ0Ptr(SP::SiconosVector newPtr)
// {
//   if(newPtr->size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setQ0Ptr: inconsistent input vector size ");
//   _q0 = newPtr;
// }

// void NewtonEulerDS::setVelocity(const SiconosVector& newValue)
// {
//   if(newValue.size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setVelocity: inconsistent input vector size ");

//   if( ! _q[1] )
//     _q[1].reset(new SimpleVector(newValue));
//   else
//     *_q[1] = newValue;
// }

// void NewtonEulerDS::setVelocityPtr(SP::SiconosVector newPtr)
// {
//   if(newPtr->size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setVelocityPtr: inconsistent input vector size ");
//   _q[1] = newPtr;
// }


// void NewtonEulerDS::setVelocity0Ptr(SP::SiconosVector newPtr)
// {
//   if(newPtr->size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setVelocity0Ptr: inconsistent input vector size ");
//   _velocity0 = newPtr;
// }

// SP::SiconosVector NewtonEulerDS::acceleration() const
// {
//   return _q[2];
// }

// void NewtonEulerDS::setP(const SiconosVector& newValue, unsigned int level)
// {
//   if(newValue.size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setP: inconsistent input vector size ");

//   if( ! _p[level] )
//     _p[level].reset(new SimpleVector(newValue));
//   else
//     *(_p[level]) = newValue;
// }

// void NewtonEulerDS::setPPtr(SP::SiconosVector newPtr, unsigned int level)
// {

//   if(newPtr->size()!= _ndof)
//     RuntimeException::selfThrow("NewtonEulerDS - setPPtr: inconsistent input vector size ");
//   _p[level] = newPtr;
// }
/*
void NewtonEulerDS::setMass(const PMMass& newValue)
{
  assert(newValue.size(0)==_ndof&&"NewtonEulerDS - setMass: inconsistent dimensions with problem size for matrix mass.");
  assert(newValue.size(1)==_ndof&&"NewtonEulerDS - setMass: inconsistent dimensions with problem size for matrix mass.");

  if( !mass  )
    mass.reset(new PMMass(newValue));
  else
    *mass = newValue;
}
void NewtonEulerDS::setFInt(const PVFint& newValue)
{
  assert(newValue.size()==_n&&"NewtonEulerDS - setFInt: inconsistent dimensions with problem size for input vector fInt");

  if( ! fInt )
    fInt.reset(new PVFint(newValue));
  else
    *_fInt = newValue;
}

void NewtonEulerDS::setFExt(const SimpleVector& newValue)
{
  assert(newValue.size()==_n&&"NewtonEulerDS - setFExt: inconsistent dimensions with problem size for input vector fExt");

  if( !fExt )
    fExt.reset(new Plugged_Vector_FTime(newValue));
  else
    *_fExt = newValue;
}

void NewtonEulerDS::setNNL(const PVNNL& newValue)
{
  assert(newValue.size()==_n&&"NewtonEulerDS - setNNL: inconsistent dimensions with problem size for input vector NNL");

  if( !NNL )
    NNL.reset(new PVNNL(newValue));
  else
    *NNL = newValue;
}

void NewtonEulerDS::setJacobianFInt(unsigned int i, const PMFint& newValue)
{
  assert(newValue.size(0)==_n&&"NewtonEulerDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");
  assert(newValue.size(1)==_n&&"NewtonEulerDS -setJacobianFInt : inconsistent dimensions with problem size for matrix JacobianFInt.");

  if( ! jacobianFInt[i] )
    jacobianFInt[i].reset(new PMFint(newValue));
  else
    *_jacobianFInt[i] = newValue;
}

void NewtonEulerDS::setJacobianNNL(unsigned int i, const PMNNL& newValue)
{
  assert(newValue.size(0)==_n&&"NewtonEulerDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");
  assert(newValue.size(1)==_n&&"NewtonEulerDS - setJacobianNNL: inconsistent dimensions with problem size for matrix JacobianNNL.");

  if( ! jacobianNNL[i] )
    jacobianNNL[i].reset(new PMNNL(newValue));
  else
    *_jacobianNNL[i] = newValue;
}
*/



// void NewtonEulerDS::computeFInt(double time){
//   if(computeFIntPtr)
//     (computeFIntPtr)(time, _n, &(*_q[0])(0), &(*_q[1])(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
// }
// void NewtonEulerDS::computeFInt(double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeFIntPtr)
//       (computeFIntPtr)(time, _n, &(*q2)(0), &(*velocity2)(0), &(*_fInt)(0), _z->size(), &(*_z)(0));
// }

void NewtonEulerDS::computeFExt(double time)
{
  if (computeFExtPtr)
    (computeFExtPtr)(time, _n, &(*_fExt)(0), _z->size(), &(*_z)(0));
}
void NewtonEulerDS::computeMmtExt(double time)
{
  if (computeMmtExtPtr)
    (computeMmtExtPtr)(time, _n, &(*_mExt)(0), _z->size(), &(*_z)(0));
}

void NewtonEulerDS::computeNNL()
{
  if (computeNNLPtr)
    (computeNNLPtr)(_n, &(*_q)(0), &(*_v)(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
}

void NewtonEulerDS::computeNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeNNLPtr)
    (computeNNLPtr)(_n, &(*q2)(0), &(*velocity2)(0), &(*_NNL)(0), _z->size(), &(*_z)(0));
  /*subtract Omega I Omega*/
}

// void NewtonEulerDS::computeJacobianqFInt(double time){
//   if(computeJacobianqFIntPtr)
//       (computeJacobianqFIntPtr)(time, _n, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianqFInt)(0,0), _z->size(), &(*_z)(0));
// }
// void NewtonEulerDS::computeJacobianqDotFInt(double time){
//   if(computeJacobianqDotFIntPtr)
//       (computeJacobianqDotFIntPtr)(time, _n, &(*_q)(0), &(*_v)(0), &(*_jacobianqDotFInt)(0,0), _z->size(), &(*_z)(0));
// }
// void NewtonEulerDS::computeJacobianZFInt(double time){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _n, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

// void NewtonEulerDS::computeJacobianqFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianqFIntPtr)
//       (computeJacobianqFIntPtr)(time, _n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianqFInt)(0,0), _z->size(), &(*_z)(0));
// }
// void NewtonEulerDS::computeJacobianqDotFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianqDotFIntPtr)
//       (computeJacobianqDotFIntPtr)(time, _n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianqDotFInt)(0,0), _z->size(), &(*_z)(0));
// }
// void NewtonEulerDS::computeJacobianZFInt( double time, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZFIntPtr)
//       (computeJacobianZFIntPtr)(time, _n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianFInt[i])(0,0), _z->size(), &(*_z)(0));
// }

void NewtonEulerDS::computeJacobianqNNL()
{
  if (computeJacobianqNNLPtr)
    (computeJacobianqNNLPtr)(_n, &(*_q)(0), &(*_v)(0), &(*_jacobianqNNL)(0, 0), _z->size(), &(*_z)(0));
}
void NewtonEulerDS::computeJacobianqDotNNL()
{
  if (computeJacobianqDotNNLPtr)
    (computeJacobianqDotNNLPtr)(_n, &(*_q)(0), &(*_v)(0), &(*_jacobianqDotNNL)(0, 0), _z->size(), &(*_z)(0));
}
// void NewtonEulerDS::computeJacobianZNNL(){
//   if(computeJacobianZNNLPtr)
//       (computeJacobianZNNLPtr)(_n, &(*_q[0])(0), &(*_q[1])(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void NewtonEulerDS::computeJacobianqNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianqNNLPtr)
    (computeJacobianqNNLPtr)(_n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianqNNL)(0, 0), _z->size(), &(*_z)(0));
}
void NewtonEulerDS::computeJacobianqDotNNL(SP::SiconosVector q2, SP::SiconosVector velocity2)
{
  if (computeJacobianqDotNNLPtr)
    (computeJacobianqDotNNLPtr)(_n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianqDotNNL)(0, 0), _z->size(), &(*_z)(0));
}
// void NewtonEulerDS::computeJacobianZNNL(unsigned int i, SP::SiconosVector q2, SP::SiconosVector velocity2){
//   if(computeJacobianZNNLPtr)
//     (computeJacobianZNNLPtr)(_n, &(*q2)(0), &(*velocity2)(0), &(*_jacobianNNL[i])(0,0), _z->size(), &(*_z)(0));
// }

void NewtonEulerDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...

  //  *_q = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_fL)
  {
    computeFL(time);
    //*_q += *_fL;
  }

  // mass and inv(mass) computatiton


  // Computes q[2] = inv(mass)*(fL+p) by solving Mq[2]=fL+p.
  // -- Case 1: if mass is constant, then a copy of imass is LU-factorized during initialization and saved into _workMatrix[invMass].
  // -- Case 2: mass is not constant, we copy it into _workMatrix[invMass]
  // Then we proceed with PLUForwardBackward.

  //  if(mass->isPlugged()) : mass may be not plugged in NewtonEulerDS children

  //  _workMatrix[invMass]->PLUForwardBackwardInPlace(*_q[2]);

}

void NewtonEulerDS::computeJacobianXRhs(double time, bool isDSup)
{
  RuntimeException::selfThrow("NewtonEulerDS::computeJacobianXRhs - not yet implemented.");
}

void NewtonEulerDS::computeFL(double time)
{
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_fL)
  {
    // 1 - Computes the required functions
    if (_fExt)
    {
      computeFExt(time);
      (boost::static_pointer_cast <SimpleVector>(_fL))->setBlock(0, *_fExt);
    }
    if (_mExt)
    {
      computeMmtExt(time);
      (boost::static_pointer_cast <SimpleVector>(_fL))->setBlock(3, *_mExt);
    }
    if (_NNL)
    {
      computeNNL();
      (*_fL) += (*_NNL);
    }
    /*computation of \Omega vectortiel I \Omega*/
    if (_I)
    {
      SimpleVector bufOmega(3);
      SimpleVector bufIOmega(3);
      SimpleVector buf(3);
      bufOmega.setValue(0, _v->getValue(3));
      bufOmega.setValue(1, _v->getValue(4));
      bufOmega.setValue(2, _v->getValue(5));
      prod(*_I, bufOmega, bufIOmega, true);
      cross_product(bufOmega, bufIOmega, buf);
      _fL->setValue(3, _fL->getValue(3) + buf.getValue(0));
      _fL->setValue(4, _fL->getValue(4) + buf.getValue(1));
      _fL->setValue(5, _fL->getValue(5) + buf.getValue(2));


    }
  }
  // else nothing.
}

// void NewtonEulerDS::computeFL(double time, SP::SiconosVector q2, SP::SiconosVector v2)
// {
//   // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
//   if( _fL )
//     {
//       //      computeFInt(time, q2, v2);
//       computeFExt(time);
//       computeMExt(time);
//       computeNNL(q2, v2);

//       // seems ok.
//       if(_fL.use_count() == 1)
//  {
//    //if not that means that fL is already (pointer-)connected with
//    // either fInt, NNL OR fExt.
//    _fL->zero();

//    if( _fInt )
//      *_fL-=*_fInt;

//    if( _fExt )
//      *_fL+=*_fExt;

//    if( _NNL )
//      *_fL-=*_NNL;
//  }
//     }
//   // else nothing.
// }

void NewtonEulerDS::computeJacobianqFL(double time)
{
  if (_jacobianqFL)
  {
    //      computeJacobianqFInt(time);
    computeJacobianqNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianqFL->zero();
      //    if( _jacobianqFInt )
      //      *_jacobianqFL-=*_jacobianqFInt;
      if (_jacobianqNNL)
        *_jacobianqFL -= *_jacobianqNNL;
    }
  }
  //else nothing.
}
void NewtonEulerDS::computeJacobianqDotFL(double time)
{
  if (_jacobianqDotFL)
  {
    //      computeJacobianqDotFInt(time);
    computeJacobianqDotNNL();

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt or jacobianNNL
      _jacobianqDotFL->zero();
      //    if( _jacobianqDotFInt )
      //      *_jacobianqDotFL-=*_jacobianqDotFInt;
      if (_jacobianqDotNNL)
        *_jacobianqDotFL -= *_jacobianqDotNNL;
    }
  }
  //else nothing.
}
// void NewtonEulerDS::computeJacobianZFL( double time){
//    RuntimeException::selfThrow("NewtonEulerDS::computeJacobianZFL - not implemented");
// }

void NewtonEulerDS::saveSpecificDataToXML()
{
  // --- other data ---
  /*  if(!dsxml)
    RuntimeException::selfThrow("NewtonEulerDS::saveDSToXML - object DynamicalSystemXML does not exist");

  SP::NewtonEulerDSXML lgptr = boost::static_pointer_cast <NewtonEulerDSXML>(dsxml);
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

      // JacobianqNNL
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

void NewtonEulerDS::display() const
{
  cout << "=====> NewtonEuler System display (number: " << _number << ")." << endl;
  cout << "- _n : " << _n << endl;
  cout << "- q " << endl;
  if (_q) _q->display();
  else cout << "-> NULL" << endl;
  cout << "- q0 " << endl;
  if (_q0) _q0->display();
  cout << "- v " << endl;
  if (_v) _v->display();
  else cout << "-> NULL" << endl;
  cout << "- v0 " << endl;
  if (_v0) _v0->display();
  else cout << "-> NULL" << endl;
  cout << "- p " << endl;
  if (_p[2]) _p[2]->display();
  else cout << "-> NULL" << endl;
  cout << "===================================== " << endl;
}

// --- Functions for memory handling ---
void NewtonEulerDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    cout << "Warning : FirstOrderNonLinearDS::initMemory with size equal to zero" << endl;
  else
  {
    _qMemory.reset(new SiconosMemory(steps));
    _vMemory.reset(new SiconosMemory(steps));
    _fLMemory.reset(new SiconosMemory(steps));
    _dotqMemory.reset(new SiconosMemory(steps));
    swapInMemory();
  }
}

void NewtonEulerDS::swapInMemory()
{

  //  _xMemory->swap(_x[0]);
  _qMemory->swap(_q);
  _vMemory->swap(_v);
  _dotqMemory->swap(_dotq);
  _fLMemory->swap(_fL);
  // initialization of the reaction force due to the non smooth law
  _p[1]->zero();
}

NewtonEulerDS* NewtonEulerDS::convert(DynamicalSystem* ds)
{
  NewtonEulerDS* lnlds = dynamic_cast<NewtonEulerDS*>(ds);
  return lnlds;
}
/*must be remove, replace by the RelativeConvergenceCriteron of the simulation*/

/*double NewtonEulerDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  SP::SiconosVector valRef = workV[NewtonSave];

  sub(*(q[0]),*valRef,*valRef);
  dsCvgIndic= valRef->norm2()/(valRef->norm2()+1);
  return (dsCvgIndic);
  }*/


void NewtonEulerDS::resetNonSmoothPart()
{
  _p[1]->zero();
}

void NewtonEulerDS::setComputeNNLFunction(const std::string& pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeNNLPtr, pluginPath, functionName);
}
void NewtonEulerDS::setComputeNNLFunction(FPtr5 fct)
{
  computeNNLPtr = fct;
}
//   void NewtonEulerDS::setComputeJacobianqFIntFunction( const std::string&  pluginPath, const std::string&  functionName){
//     Plugin::setFunction(&computeJacobianqFIntPtr, pluginPath,functionName);
//   }
//   void NewtonEulerDS::setComputeJacobianqDotFIntFunction( const std::string&  pluginPath, const std::string&  functionName){
//     Plugin::setFunction(&computeJacobianqDotFIntPtr, pluginPath,functionName);
//   }
//   void NewtonEulerDS::setComputeJacobianqFIntFunction(FPtr6 fct){computeJacobianqFIntPtr=fct;}
//   void NewtonEulerDS::setComputeJacobianqDotFIntFunction(FPtr6 fct){computeJacobianqDotFIntPtr=fct;}
void NewtonEulerDS::setComputeJacobianqNNLFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianqNNLPtr, pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianqDotNNLFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  Plugin::setFunction(&computeJacobianqDotNNLPtr, pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianqNNLFunction(FPtr5 fct)
{
  computeJacobianqNNLPtr = fct;
}
void NewtonEulerDS::setComputeJacobianqDotNNLFunction(FPtr5 fct)
{
  computeJacobianqDotNNLPtr = fct;
}

void NewtonEulerDS::updateT()
{
  double q0 = _q->getValue(3) / 2.0;
  double q1 = _q->getValue(4) / 2.0;
  double q2 = _q->getValue(5) / 2.0;
  double q3 = _q->getValue(6) / 2.0;
  _T->setValue(3, 3, -q1);
  _T->setValue(3, 4, -q2);
  _T->setValue(3, 5, -q3);
  _T->setValue(4, 3, q0);
  _T->setValue(4, 4, -q3);
  _T->setValue(4, 5, q2);
  _T->setValue(5, 3, q3);
  _T->setValue(5, 4, q0);
  _T->setValue(5, 5, -q1);
  _T->setValue(6, 3, -q2);
  _T->setValue(6, 4, q1);
  _T->setValue(6, 5, q0);

}
