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
#include "NewtonEulerDS.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"
#include <boost/math/quaternion.hpp>

#include <iostream>
// #define DEBUG_NOCOLOR
// //#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include <debug.h>


void computeRotationMatrix(double q0, double q1, double q2, double q3,
                           SP::SimpleMatrix rotationMatrix)
{

  /* Brute force version by multiplication of quaternion
   */
  // ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  // ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  // ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  // ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  // ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  // ::boost::math::quaternion<double>    quatBuff;
  // quatBuff = quatQ * quatx * quatcQ;
  // rotationMatrix->setValue(0, 0, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 0, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 0, quatBuff.R_component_4());
  // quatBuff = quatQ * quaty * quatcQ;
  // rotationMatrix->setValue(0, 1, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 1, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 1, quatBuff.R_component_4());
  // quatBuff = quatQ * quatz * quatcQ;
  // rotationMatrix->setValue(0, 2, quatBuff.R_component_2());
  // rotationMatrix->setValue(1, 2, quatBuff.R_component_3());
  // rotationMatrix->setValue(2, 2, quatBuff.R_component_4());

  /* direct computation https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation */
  rotationMatrix->setValue(0, 0,     q0*q0 +q1*q1 -q2*q2 -q3*q3 );
  rotationMatrix->setValue(0, 1, 2.0*(q1*q2        - q0*q3));
  rotationMatrix->setValue(0, 2, 2.0*(q1*q3        + q0*q2));

  rotationMatrix->setValue(1, 0, 2.0*(q1*q2        + q0*q3) );
  rotationMatrix->setValue(1, 1,     q0*q0 -q1*q1 +q2*q2 -q3*q3 );
  rotationMatrix->setValue(1, 2, 2.0*(q2*q3        - q0*q1) );

  rotationMatrix->setValue(2, 0, 2.0*(q1*q3        - q0*q2) );
  rotationMatrix->setValue(2, 1, 2.0*(q2*q3         + q0*q1) );
  rotationMatrix->setValue(2, 2,     q0*q0 -q1*q1 -q2*q2 +q3*q3 );
}

void rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SiconosVector v )
{
  DEBUG_BEGIN("::rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SiconosVector v )\n");
  DEBUG_EXPR(v->display(););
  DEBUG_PRINTF("( q0 = %e,  q1 = %e,  q2= %e,  q3= %e )\n", q0,q1,q2,q3);

  // First way. Using the rotation matrix
  // SP::SimpleMatrix rotationMatrix(new SimpleMatrix(3,3));
  // SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,q1,q2,q3, rotationMatrix);
  // prod(*rotationMatrix, *v,  tmp);
  // *v =tmp;

  // Second way. Using the transpose of the rotation matrix
  // SP::SimpleMatrix rotationMatrix(new SimpleMatrix(3,3));
  // SiconosVector tmp(3);
  // ::computeRotationMatrix(q0,-q1,-q2,-q3, rotationMatrix);
  // prod(*v, *rotationMatrix,  tmp);
  // *v =tmp;

  // Third way. cross product and axis angle
  // see http://www.geometrictools.com/Documentation/RotationIssues.pdf
  // SP::SiconosVector axis(new SiconosVector(3));
  // double angle = ::getAxisAngle(q0,q1,q2,q3, axis);
  // SiconosVector t(3), tmp(3);
  // cross_product(*axis,*v,t);
  // cross_product(*axis,t,tmp);
  // *v += sin(angle)*t + (1.0-cos(angle))*tmp;

  // Direct computation with cross product
  SiconosVector t(3), tmp(3);
  SiconosVector qvect(3);
  qvect(0)=q1;
  qvect(1)=q2;
  qvect(2)=q3;
  cross_product(qvect,*v,t);
  t *= 2.0;
  cross_product(qvect,t,tmp);
  *v += tmp;
  *v += q0*t;
  DEBUG_EXPR(v->display(););
  DEBUG_END("::rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SiconosVector v )\n");
}

void rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SimpleMatrix m )
{
  DEBUG_BEGIN("::rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SimpleMatrix m )\n");
  DEBUG_EXPR(m->display(););
  DEBUG_PRINTF("( q0 = %e,  q1 = %e,  q2= %e,  q3= %e )\n", q0,q1,q2,q3);

  // Direct computation with cross product for each column
  assert(m->size(0) == 3 && "::rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SimpleMatrix m ) m must have 3 rows" );
  SiconosVector v(3);
  SiconosVector t(3), tmp(3);
  SiconosVector qvect(3);
  qvect(0)=q1;
  qvect(1)=q2;
  qvect(2)=q3;
  for(unsigned int j = 0; j < m->size(1); j++)
  {
    v(0) = m->getValue(0,j);
    v(1) = m->getValue(1,j);
    v(2) = m->getValue(2,j);
    cross_product(qvect,v,t);
    t *= 2.0;
    cross_product(qvect,t,tmp);
    v += tmp;
    v += q0*t;
    m->setValue(0,j,v(0));
    m->setValue(1,j,v(1));
    m->setValue(2,j,v(2));
  }
  DEBUG_EXPR(m->display(););
  DEBUG_END("::rotateAbsToBody(double q0, double q1, double q2, double q3, SP::SimpleMatrix m )\n");
}


void rotateAbsToBody(SP::SiconosVector q, SP::SiconosVector v )
{
  DEBUG_BEGIN("::rotateAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::rotateAbsToBody(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), v);
  DEBUG_END("::rotateAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
}

void rotateAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )
{
  DEBUG_BEGIN("::rotateAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::rotateAbsToBody(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),m);
  DEBUG_END("::rotateAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m)\n");
}

void changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v )
{
  DEBUG_BEGIN("::changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::rotateAbsToBody(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6), v);
  DEBUG_END("::changeFrameAbsToBody(SP::SiconosVector q, SP::SiconosVector v )\n");
}
void changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )
{
  DEBUG_BEGIN("::changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::rotateAbsToBody(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6), m);
  DEBUG_END("::changeFrameAbsToBody(SP::SiconosVector q, SP::SimpleMatrix m )\n");
}

void changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v )
{
  DEBUG_BEGIN("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v )\n");
  ::rotateAbsToBody(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), v);
  DEBUG_END("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SiconosVector v )\n");
}
void changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m )
{
  DEBUG_BEGIN("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m )\n");
  ::rotateAbsToBody(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6), m);
  DEBUG_END("::changeFrameBodyToAbs(SP::SiconosVector q, SP::SimpleMatrix m )\n");
}



void computeRotationMatrix(SP::SiconosVector q, SP::SimpleMatrix rotationMatrix  )
{
  ::computeRotationMatrix(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),
                          rotationMatrix);
}
void computeRotationMatrixTransposed(SP::SiconosVector q, SP::SimpleMatrix rotationMatrix )
{
  ::computeRotationMatrix(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6),
                               rotationMatrix);
}

double getAxisAngle(double q0, double q1, double q2, double q3, SP::SiconosVector axis )
{
  DEBUG_BEGIN("getAxisAngle(double q0, double q1, double q2, double q3, SP::SiconosVector axis )\n");
  double angle = acos(q0) *2.0;
  //double f = sin( angle *0.5);
  double f = sqrt(1-q0*q0); // cheaper than sin ?
  if (f !=0.0)
  {
    axis->setValue(0, q1/f);
    axis->setValue(1, q2/f);
    axis->setValue(2, q3/f);
  }
  else
  {
    axis->zero();
  }
  DEBUG_PRINTF("angle= %12.8e\n", angle);
  DEBUG_EXPR(axis->display(););
  DEBUG_END("getAxisAngle(double q0, double q1, double q2, double q3, SP::SiconosVector axis )\n");
  return angle;
}
double getAxisAngle(SP::SiconosVector q, SP::SiconosVector axis )
{
  double angle = ::getAxisAngle(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6),axis );
  return angle;
}


void setAxisAngle(SP::SiconosVector q, SP::SiconosVector axis, double angle  )
{

  double normv = sqrt(axis->getValue(0)*axis->getValue(0)+
                      axis->getValue(1)*axis->getValue(1)+
                      axis->getValue(2)*axis->getValue(2));

  if (normv != 0.0)
  {
    q->setValue(3, cos(angle/2.0));
    double f = sin( angle *0.5);
    q->setValue(4, axis->getValue(0)/normv * f);
    q->setValue(5, axis->getValue(1)/normv * f);
    q->setValue(6, axis->getValue(2)/normv * f);
  }
  else
  {
    // we assume that the axis has to be well defined.
    // we return a non unit quaternion.
    q->zero();
  }
}

void normalizeq(SP::SiconosVector q)
{
  double normq = sqrt(q->getValue(3) * q->getValue(3) +
                      q->getValue(4) * q->getValue(4) +
                      q->getValue(5) * q->getValue(5) +
                      q->getValue(6) * q->getValue(6));
  assert(normq > 0);
  normq = 1.0 / normq;
  q->setValue(3, q->getValue(3) * normq);
  q->setValue(4, q->getValue(4) * normq);
  q->setValue(5, q->getValue(5) * normq);
  q->setValue(6, q->getValue(6) * normq);
}


void computeT(SP::SiconosVector q, SP::SimpleMatrix T)
{
  DEBUG_BEGIN("computeT(SP::SiconosVector q, SP::SimpleMatrix T)\n")
  //  std::cout <<"\n NewtonEulerDS::computeT(SP::SiconosVector q)\n  " <<std::endl;
  double q0 = q->getValue(3) / 2.0;
  double q1 = q->getValue(4) / 2.0;
  double q2 = q->getValue(5) / 2.0;
  double q3 = q->getValue(6) / 2.0;
  T->setValue(3, 3, -q1);
  T->setValue(3, 4, -q2);
  T->setValue(3, 5, -q3);
  T->setValue(4, 3, q0);
  T->setValue(4, 4, -q3);
  T->setValue(4, 5, q2);
  T->setValue(5, 3, q3);
  T->setValue(5, 4, q0);
  T->setValue(5, 5, -q1);
  T->setValue(6, 3, -q2);
  T->setValue(6, 4, q1);
  T->setValue(6, 5, q0);
  DEBUG_END("computeT(SP::SiconosVector q, SP::SimpleMatrix T)\n")

}



// Private function to set linked with members of Dynamical top class
void NewtonEulerDS::connectToDS()
{
  // dim
  _n = 2 * 3;

}




// From a set of data; Mass filled-in directly from a siconosMatrix -
// This constructor leads to the minimum NewtonEuler System form: \f$ M\ddot q = p \f$
/*
Q0 : contains the center of mass coordinate, and the quaternion initial. (dim(Q0)=7)
Velocity0 : contains the initial velocity of center of mass and the omega initial. (dim(Velocity0)=6)
*/
NewtonEulerDS::NewtonEulerDS(): DynamicalSystem(6),
                                nullifyFGyr(false),
                                _computeJacobianFIntqByFD(false),
                                _computeJacobianFIntvByFD(false),
                                _computeJacobianMIntqByFD(false),
                                _computeJacobianMIntvByFD(false),
                                _epsilonFD(sqrt(std::numeric_limits< double >::epsilon()))
{
  _p.resize(3);
  _p[0].reset(new SiconosVector());
  _p[1].reset(new SiconosVector(_n)); // Needed in NewtonEulerR
  _p[2].reset(new SiconosVector());
  zeroPlugin();
  //assert(0);

  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  _qDim = 7;
  _n = 6;

  // Current state
  _q.reset(new SiconosVector(_qDim));
  // _deltaq.reset(new SiconosVector(_qDim));
  _v.reset(new SiconosVector(_n));

  _dotq.reset(new SiconosVector(_qDim));
  _workspace[freeresidu].reset(new SiconosVector(_n));
  _workspace[free].reset(new SiconosVector(dimension()));
  _massMatrix.reset(new SimpleMatrix(_n, _n));
  _luW.reset(new SimpleMatrix(_n, _n));
  _massMatrix->zero();
  _T.reset(new SimpleMatrix(_qDim, _n));

  _scalarMass = 0.;
}

void NewtonEulerDS::internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0,
                                 double mass , SP::SiconosMatrix inertialMatrix)
{
  DEBUG_BEGIN("NewtonEulerDS::internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix)\n");
  _p.resize(3);
  _p[0].reset(new SiconosVector());
  _p[1].reset(new SiconosVector(_n)); // Needed in NewtonEulerR
  _p[2].reset(new SiconosVector());
  zeroPlugin();
  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  _scalarMass = mass;
  _qDim = 7;
  _n = 6;

  // Initial conditions
  _q0 = Q0;
  _v0 = Velocity0;


  // Current state
  _q.reset(new SiconosVector(_qDim));
  // _deltaq.reset(new SiconosVector(_qDim));
  _v.reset(new SiconosVector(_n));
  (*_q) = (*_q0);
  _dotq.reset(new SiconosVector(_qDim));
  _massMatrix.reset(new SimpleMatrix(_n, _n));
  _jacobianFGyrv.reset(new SimpleMatrix(_n, _n));
  _luW.reset(new SimpleMatrix(_n, _n));
  if (inertialMatrix)
    _I = inertialMatrix;
  else {
    _I.reset(new SimpleMatrix(3, 3));
    _I->zero();
    _I->setValue(0,0, 1.0);
    _I->setValue(1,1, 1.0);
    _I->setValue(2,2, 1.0);
  }
  updateMassMatrix();
  _workspace[freeresidu].reset(new SiconosVector(_n));
  _workspace[free].reset(new SiconosVector(dimension()));

  _T.reset(new SimpleMatrix(_qDim, _n));
  _T->zero();
  _T->setValue(0, 0, 1.0);
  _T->setValue(1, 1, 1.0);
  _T->setValue(2, 2, 1.0);
  computeT();
  initForces();
  DEBUG_END("NewtonEulerDS::internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix)\n");
}

void NewtonEulerDS::updateMassMatrix()
{
  _massMatrix->zero();
  _massMatrix->setValue(0, 0, _scalarMass);
  _massMatrix->setValue(1, 1, _scalarMass);
  _massMatrix->setValue(2, 2, _scalarMass);
  Index dimIndex(2);
  dimIndex[0] = 3;
  dimIndex[1] = 3;
  Index startIndex(4);
  startIndex[0] = 0;
  startIndex[1] = 0;
  startIndex[2] = 3;
  startIndex[3] = 3;
  setBlock(_I, _massMatrix, dimIndex, startIndex);
}

NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Velocity0,
                             double  mass, SP::SiconosMatrix inertialMatrix):
  DynamicalSystem(6),
  _computeJacobianFIntqByFD(false),
  _computeJacobianFIntvByFD(false),
  _computeJacobianMIntqByFD(false),
  _computeJacobianMIntvByFD(false),
  _epsilonFD(sqrt(std::numeric_limits< double >::epsilon()))
{
  internalInit(Q0, Velocity0, mass, inertialMatrix);
}

void NewtonEulerDS::zeroPlugin()
{
  _pluginFExt.reset(new PluggedObject());
  _pluginMExt.reset(new PluggedObject());
  _pluginFInt.reset(new PluggedObject());
  _pluginMInt.reset(new PluggedObject());
  _pluginJacqFInt.reset(new PluggedObject());
  _pluginJacvFInt.reset(new PluggedObject());
  _pluginJacqMInt.reset(new PluggedObject());
  _pluginJacvMInt.reset(new PluggedObject());
}

// Destructor
NewtonEulerDS::~NewtonEulerDS()
{
}

void NewtonEulerDS::setInertia(double ix, double iy, double iz)
{
  _I->zero();

  (*_I)(0, 0) = ix;
  (*_I)(1, 1) = iy;
  (*_I)(2, 2) = iz;

  updateMassMatrix();
}

bool NewtonEulerDS::checkDynamicalSystem()
{
  bool output = true;
  // ndof

  // q0 and velocity0
  if (! _q0 || ! _v0)
  {
    RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }


  // fInt
  //   if( ( _fInt && computeFIntPtr) && ( ! _jacobianFIntq || ! _jacobianFIntv ) )
  //     // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  //     {
  //       RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
  //       output = false;
  //     }


  if (!output) std::cout << "NewtonEulerDS Warning: your dynamical system seems to be uncomplete (check = false)" <<std::endl;
  return output;
}

void NewtonEulerDS::initializeNonSmoothInput(unsigned int level)
{
  DEBUG_PRINTF("NewtonEulerDS::initializeNonSmoothInput(unsigned int level) for level = %i\n",level);

  if (_p[level]->size() == 0)
  {
    if (level == 0)
    {
      _p[0]->resize(_qDim);
    }
    else
    {
      _p[level]->resize(_n);
    }
  }


#ifdef DEBUG_MESSAGES
  DEBUG_PRINT("display() after initialization");
  display();
#endif
}

void NewtonEulerDS::initForces()
{
  DEBUG_BEGIN("NewtonEulerDS::initForces()\n")
  _forces.reset(new SiconosVector(_n));
  _fGyr.reset(new SiconosVector(3,0.0));
  _jacobianFGyrv.reset(new SimpleMatrix(_n, _n));
  _jacobianvForces.reset(new SimpleMatrix(_n, _n));
  DEBUG_END("NewtonEulerDS::initForces()\n")
}

void NewtonEulerDS::initRhs(double time)
{
  //  _workMatrix.resize(sizeWorkMat);

  // Solve Mq[2]=fL+p.
  //*_q = *(_p[2]); // Warning: r/p update is done in Interactions/Relations

  if (_forces)
  {
    computeForces(time);
    //      *_q += *_forces;
  }

}

void NewtonEulerDS::initialize(double time, unsigned int sizeOfMemory)
{
  // set q and q[1] to q0 and velocity0, initialize acceleration.
  *_q = *_q0;
  *_v = *_v0;

  // If z has not been set, we initialize it with a null vector of size 1, since z is required in plug-in functions call.
  if (! _z)
    _z.reset(new SiconosVector(1));

  if (_pluginFExt->fPtr && !_fExt)
    _fExt.reset(new SiconosVector(3, 0));

  if (_pluginMExt->fPtr && !_mExt)
  {
    _mExt.reset(new SiconosVector(3, 0));
  }
  if (_pluginFInt->fPtr && !_fInt)
    _fInt.reset(new SiconosVector(3, 0));

  if ((_pluginJacqFInt->fPtr  || _computeJacobianFIntqByFD) && !_jacobianFIntq)
  {
    _jacobianFIntq.reset(new SimpleMatrix(3, _qDim));
    if (!_jacobianqForces)
      _jacobianqForces.reset(new SimpleMatrix(_n, _qDim));
  }

  if ((_pluginJacvFInt->fPtr || _computeJacobianFIntvByFD) && !_jacobianFIntv)
    _jacobianFIntv.reset(new SimpleMatrix(3, _n));

  if (_pluginMInt->fPtr && !_mInt)
    _mInt.reset(new SiconosVector(3, 0));

  if ((_pluginJacqMInt->fPtr || _computeJacobianMIntqByFD) && !_jacobianMIntq)
  {
    if (!_jacobianqForces)
      _jacobianqForces.reset(new SimpleMatrix(_n, _qDim));
    _jacobianMIntq.reset(new SimpleMatrix(3, _qDim));
  }
  if ((_pluginJacvMInt->fPtr || _computeJacobianMIntvByFD) && !_jacobianMIntv)
    _jacobianMIntv.reset(new SimpleMatrix(3, _n));

  if(_mExt)
  {
    _mExtObj.reset(new SiconosVector(3, 0));
     if (!_jacobianqForces)
       _jacobianqForces.reset(new SimpleMatrix(_n, _qDim));
     _jacobianMExtObjq.reset(new SimpleMatrix(3, _qDim));
  }

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links.
  // No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  checkDynamicalSystem();

  initRhs(time);


  if (_boundaryConditions)
  {
    _reactionToBoundaryConditions.reset(new SiconosVector(_boundaryConditions->velocityIndices()->size()));
  }

  // Initialize memory vectors
  initMemory(sizeOfMemory);

}

void NewtonEulerDS::computeFExt(double time)
{
  if (_pluginFExt->fPtr)
    ((FExt_NE)_pluginFExt->fPtr)(time, &(*_fExt)(0), _qDim, &(*_q0)(0) ); // parameter z are assumed to be equal to q0
}

void NewtonEulerDS::computeMExt(double time)
{
  computeMExt(time,  _mExt);
}
void NewtonEulerDS::computeMExt(double time, SP::SiconosVector mExt)
{
  // Take the constant value of exists otherwise call the plugin
  *mExt = *_mExt;
  if (_pluginMExt->fPtr)
    ((FExt_NE)_pluginMExt->fPtr)(time, &(*mExt)(0), _qDim, &(*_q0)(0) ); // parameter z are assumed to be equal to q0
}
void NewtonEulerDS::computeMExtObj(double time)
{
  computeMExtObj(time, _q,  _mExtObj );
}

void NewtonEulerDS::computeMExtObj(double time, SP::SiconosVector q,  SP::SiconosVector mExtObj )
{
  DEBUG_BEGIN("NewtonEulerDS::computeMExtObj(...)\n");

   // Call plugin if needed
  computeMExt(time, mExtObj); // Compute Mext and store it direclty in MExtObj
  ::changeFrameAbsToBody(q,mExtObj);

  DEBUG_EXPR(mExtObj->display());
  DEBUG_END("NewtonEulerDS::computeMExtObj(...)\n");
}





void NewtonEulerDS::computeJacobianMExtObjqByFD(double time, SP::SiconosVector q)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianMExtvByFD(...)\n");
  SP::SiconosVector mExtObj(new SiconosVector(3));
  computeMExtObj(time, q, mExtObj);

  double mExtObj0 = mExtObj->getValue(0);
  double mExtObj1 = mExtObj->getValue(1);
  double mExtObj2 = mExtObj->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));
  _jacobianMExtObjq->zero();
  (*qeps)(0) += _epsilonFD;
  for (int j =0; j < 7; j++)
  {
    computeMExtObj(time, qeps, mExtObj);
    _jacobianMExtObjq->setValue(0,j,  (mExtObj->getValue(0) - mExtObj0)/_epsilonFD );
    _jacobianMExtObjq->setValue(1,j,  (mExtObj->getValue(1) - mExtObj1)/_epsilonFD );
    _jacobianMExtObjq->setValue(2,j,  (mExtObj->getValue(2) - mExtObj2)/_epsilonFD );
    (*qeps)(j) -= _epsilonFD;
    if (j<6) (*qeps)(j+1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianMExtObjq->display(););
  DEBUG_END("NewtonEulerDS::computeJacobianMExtvByFD(...)\n");


}

void NewtonEulerDS::computeFInt(double time)
{
  computeFInt(time, _q, _v);
}


void NewtonEulerDS::computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v)
{
  computeFInt(time,  q,  v, _fInt);
}

void NewtonEulerDS::computeFInt(double time, SP::SiconosVector q, SP::SiconosVector v, SP::SiconosVector fInt)
{
  if (_pluginFInt->fPtr)
    ((FInt_NE)_pluginFInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*fInt)(0), _qDim,  &(*_q0)(0));// parameter z are assumed to be equal to q0
}


void NewtonEulerDS::computeMInt(double time)
{
  computeMInt(time, _q, _v);
}


void NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v)
{
  computeMInt(time, q, v, _mInt);
}

void NewtonEulerDS::computeMInt(double time, SP::SiconosVector q, SP::SiconosVector v, SP::SiconosVector mInt)
{
  if (_pluginMInt->fPtr)
    ((FInt_NE)_pluginMInt->fPtr)(time, &(*q)(0), &(*v)(0), &(*mInt)(0), _qDim,  &(*_q0)(0));// parameter z are assumed to be equal to q0
}

void NewtonEulerDS::computeJacobianFIntq(double time)
{
  computeJacobianFIntq(time, _q, _v);
}
void NewtonEulerDS::computeJacobianFIntv(double time)
{
  computeJacobianFIntv(time, _q, _v);
}

void NewtonEulerDS::computeJacobianFIntq(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_PRINT("NewtonEulerDS::computeJacobianFIntq(...) starts");
  if (_pluginJacqFInt->fPtr)
    ((FInt_NE)_pluginJacqFInt->fPtr)(time, &(*q)(0), &(*velocity)(0), &(*_jacobianFIntq)(0, 0), _qDim,  &(*_q0)(0));
  else if (_computeJacobianFIntqByFD)
    computeJacobianFIntqByFD(time, q, velocity);
  DEBUG_EXPR(_jacobianFIntq->display(););
  DEBUG_END("NewtonEulerDS::computeJacobianFIntq(...)");
}

void NewtonEulerDS::computeJacobianFIntqByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFIntqByFD(...)\n");
  SP::SiconosVector fInt(new SiconosVector(3));
  computeFInt(time, q, velocity, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));
  _jacobianFIntq->zero();
  (*qeps)(0) += _epsilonFD;
  for (int j =0; j < 7; j++)
  {
    computeFInt(time, qeps, velocity, fInt);
    _jacobianFIntq->setValue(0,j,  (fInt->getValue(0) - fInt0)/_epsilonFD );
    _jacobianFIntq->setValue(1,j,  (fInt->getValue(1) - fInt1)/_epsilonFD );
    _jacobianFIntq->setValue(2,j,  (fInt->getValue(2) - fInt2)/_epsilonFD );
    (*qeps)(j) -= _epsilonFD;
    if (j<6) (*qeps)(j+1) += _epsilonFD;
  }
  DEBUG_END("NewtonEulerDS::computeJacobianFIntqByFD(...)\n");


}

void NewtonEulerDS::computeJacobianFIntv(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  if (_pluginJacvFInt->fPtr)
    ((FInt_NE)_pluginJacvFInt->fPtr)(time, &(*q)(0), &(*velocity)(0), &(*_jacobianFIntv)(0, 0), _qDim,  &(*_q0)(0));
  else if (_computeJacobianFIntvByFD)
    computeJacobianFIntvByFD(time, q, velocity);
}

void NewtonEulerDS::computeJacobianFIntvByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFIntvByFD(...)\n");
  SP::SiconosVector fInt(new SiconosVector(3));
  computeFInt(time, q, velocity, fInt);

  double fInt0 = fInt->getValue(0);
  double fInt1 = fInt->getValue(1);
  double fInt2 = fInt->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*velocity));
  _jacobianFIntv->zero();

  (*veps)(0) += _epsilonFD;
  for (int j =0; j < 6; j++)
  {
    computeFInt(time, q, veps, fInt);
    _jacobianFIntv->setValue(0,j,  (fInt->getValue(0) - fInt0)/_epsilonFD );
    _jacobianFIntv->setValue(1,j,  (fInt->getValue(1) - fInt1)/_epsilonFD );
    _jacobianFIntv->setValue(2,j,  (fInt->getValue(2) - fInt2)/_epsilonFD );
    (*veps)(j) -= _epsilonFD;
    if (j<5) (*veps)(j+1) += _epsilonFD;
  }

  DEBUG_END("NewtonEulerDS::computeJacobianFIntvByFD(...)\n");


}
void NewtonEulerDS::computeJacobianFGyrvByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFGyrvByFD(...)\n");
  SP::SiconosVector fGyr(new SiconosVector(3));
  computeFGyr(velocity, fGyr);

  double fGyr0 = fGyr->getValue(0);
  double fGyr1 = fGyr->getValue(1);
  double fGyr2 = fGyr->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*velocity));
  _jacobianFGyrv->zero();


  (*veps)(0) += _epsilonFD;
  for (int j =0; j < 6; j++)
  {
    computeFGyr(veps, fGyr);
    _jacobianFGyrv->setValue(3,j,  (fGyr->getValue(0) - fGyr0)/_epsilonFD );
    _jacobianFGyrv->setValue(4,j,  (fGyr->getValue(1) - fGyr1)/_epsilonFD );
    _jacobianFGyrv->setValue(5,j,  (fGyr->getValue(2) - fGyr2)/_epsilonFD );
    (*veps)(j) -= _epsilonFD;
    if (j<5) (*veps)(j+1) += _epsilonFD;
  }
  DEBUG_EXPR(_jacobianFGyrv->display());
  DEBUG_END("NewtonEulerDS::computeJacobianFGyrvByFD(...)\n");


}
void NewtonEulerDS::computeJacobianMIntq(double time)
{
  computeJacobianMIntq(time, _q, _v);
}
void NewtonEulerDS::computeJacobianMIntv(double time)
{
  computeJacobianMIntv(time, _q, _v);
}

void NewtonEulerDS::computeJacobianMIntq(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntq(...) starts");
  if (_pluginJacqMInt->fPtr)
    ((FInt_NE)_pluginJacqMInt->fPtr)(time, &(*q)(0), &(*velocity)(0), &(*_jacobianMIntq)(0, 0), _qDim,  &(*_q0)(0));
  else if (_computeJacobianMIntqByFD)
    computeJacobianMIntqByFD(time, q, velocity);
  DEBUG_EXPR(_jacobianMIntq->display());
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntq(...) ends");

}

void NewtonEulerDS::computeJacobianMIntqByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntqByFD(...) starts\n");

  SP::SiconosVector mInt(new SiconosVector(3));
  computeMInt(time, q, velocity, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  SP::SiconosVector qeps(new SiconosVector(*q));

  (*qeps)(0) += _epsilonFD;
  for (int j =0; j < 7; j++)
  {
    computeMInt(time, qeps, velocity, mInt);
    _jacobianMIntq->setValue(0,j,  (mInt->getValue(0) - mInt0)/_epsilonFD );
    _jacobianMIntq->setValue(1,j,  (mInt->getValue(1) - mInt1)/_epsilonFD );
    _jacobianMIntq->setValue(2,j,  (mInt->getValue(2) - mInt2)/_epsilonFD );
    (*qeps)(j) -= _epsilonFD;
    if (j<6) (*qeps)(j+1) += _epsilonFD;
  }
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntqByFD(...) ends\n");
}

void NewtonEulerDS::computeJacobianMIntv(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  if (_pluginJacvMInt->fPtr)
    ((FInt_NE)_pluginJacvMInt->fPtr)(time, &(*q)(0), &(*velocity)(0), &(*_jacobianMIntv)(0, 0), _qDim,  &(*_q0)(0));
  else if (_computeJacobianMIntvByFD)
    computeJacobianMIntvByFD(time,  q, velocity);
}

void NewtonEulerDS::computeJacobianMIntvByFD(double time, SP::SiconosVector q, SP::SiconosVector velocity)
{
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntvByFD(...) starts\n");

  SP::SiconosVector mInt(new SiconosVector(3));
  computeMInt(time, q, velocity, mInt);
  double mInt0 = mInt->getValue(0);
  double mInt1 = mInt->getValue(1);
  double mInt2 = mInt->getValue(2);

  SP::SiconosVector veps(new SiconosVector(*velocity));

  (*veps)(0) += _epsilonFD;
  for (int j =0; j < 6; j++)
  {
    computeMInt(time, q, veps, mInt);
    _jacobianMIntv->setValue(0,j,  (mInt->getValue(0) - mInt0)/_epsilonFD );
    _jacobianMIntv->setValue(1,j,  (mInt->getValue(1) - mInt1)/_epsilonFD );
    _jacobianMIntv->setValue(2,j,  (mInt->getValue(2) - mInt2)/_epsilonFD );
    (*veps)(j) -= _epsilonFD;
    if (j<5) (*veps)(j+1) += _epsilonFD;
  }
  DEBUG_PRINT("NewtonEulerDS::computeJacobianMIntvByFD(...) ends\n");
}


void NewtonEulerDS::computeRhs(double time, bool isDSup)
{
  // if isDSup == true, this means that there is no need to re-compute mass ...
  //  *_q = *(_p[2]); // Warning: r/p update is done in Interactions/Relations
  if (_forces)
  {
    computeForces(time);
    //*_q += *_forces;
  }
}

void NewtonEulerDS::computeJacobianRhsx(double time, bool isDSup)
{
  RuntimeException::selfThrow("NewtonEulerDS::computeJacobianRhsx - not yet implemented.");
}

void NewtonEulerDS::computeForces(double time)
{
  computeForces(time, _q, _v);
}

void NewtonEulerDS::computeFGyr(SP::SiconosVector v, SP::SiconosVector fGyr)
{
  // computation of \Omega times I \Omega (FGyr is in the l.h.s of the equation of motion)
  DEBUG_BEGIN("NewtonEulerDS::computeFGyr(SP::SiconosVector v, SP::SiconosVector fGyr)\n");
  if (_I)
  {
    DEBUG_EXPR( _I->display());
    DEBUG_EXPR( v->display());
    SiconosVector omega(3);
    SiconosVector iomega(3);
    omega.setValue(0, v->getValue(3));
    omega.setValue(1, v->getValue(4));
    omega.setValue(2, v->getValue(5));
    prod(*_I, omega, iomega, true);
    cross_product(omega, iomega, *fGyr);
  }
  DEBUG_EXPR(fGyr->display());
  DEBUG_END("NewtonEulerDS::computeFGyr(SP::SiconosVector v, SP::SiconosVector fGyr)\n");

}
void NewtonEulerDS::computeFGyr(SP::SiconosVector v)
{
  /*computation of \Omega times I \Omega*/
  //DEBUG_BEGIN("NewtonEulerDS::computeFGyr(SP::SiconosVector v)\n");
  computeFGyr( v, _fGyr);
  //DEBUG_END("NewtonEulerDS::computeFGyr(SP::SiconosVector v)\n");

}


void NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector v)
{
  DEBUG_BEGIN("NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector v)\n")
  // Warning: an operator (fInt ...) may be set (ie allocated and not NULL) but not plugged, that's why two steps are required here.
  if (_forces)
  {
    _forces->zero();
    // 1 - Computes the required functions
    if (_fExt)
    {
      computeFExt(time);
      _forces->setBlock(0, *_fExt);
    }
    if (_mExt)
    {
      computeMExtObj(time,q,_mExtObj);
      _forces->setBlock(3, *_mExtObj);
    }
    if (_fInt)
    {
      computeFInt(time, q, v);
      // std::cout << "_fInt : "<< std::endl;
      // _fInt->display();
      _forces->setValue(0, _forces->getValue(0) - _fInt->getValue(0));
      _forces->setValue(1, _forces->getValue(1) - _fInt->getValue(1));
      _forces->setValue(2, _forces->getValue(2) - _fInt->getValue(2));

    }
    if (_mInt)
    {
      computeMInt(time, q , v);
      SiconosVector aux(3);
      ::changeFrameAbsToBody(q,_mInt); // We have to be sure that MInt is expressed in inertial frame.
      _forces->setValue(3, _forces->getValue(3) - _mInt->getValue(0));
      _forces->setValue(4, _forces->getValue(4) - _mInt->getValue(1));
      _forces->setValue(5, _forces->getValue(5) - _mInt->getValue(2));
    }
    if (!nullifyFGyr)
    {
      computeFGyr(v);
      _forces->setValue(3, _forces->getValue(3) - _fGyr->getValue(0));
      _forces->setValue(4, _forces->getValue(4) - _fGyr->getValue(1));
      _forces->setValue(5, _forces->getValue(5) - _fGyr->getValue(2));
    }
    DEBUG_EXPR(_forces->display());
    DEBUG_END("NewtonEulerDS::computeForces(double time, SP::SiconosVector q, SP::SiconosVector v)\n")

    // std::cout << "_forces : "<< std::endl;
    // _forces->display();
  }
  // else nothing.
}

void NewtonEulerDS::computeJacobianqForces(double time)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianqForces(double time) \n");
  if (_jacobianqForces)
  {
    _jacobianqForces->zero();
    if (_jacobianFIntq)
    {
      computeJacobianFIntq(time);
      _jacobianqForces->setBlock(0,0,-1.0 * *_jacobianFIntq);
    }
    if (_jacobianMIntq)
    {
      computeJacobianMIntq(time);
      ::changeFrameAbsToBody(_q,_jacobianMIntq); // We have to ensure that the Jacobian is expressed in inertial frame
    }
    if (_jacobianMExtObjq)
    {
      computeJacobianMExtObjqByFD(time, _q);
      _jacobianqForces->setBlock(3,0,1.0* *_jacobianMExtObjq);
    }
    DEBUG_EXPR(_jacobianqForces->display(););
  }
  //else nothing.
  DEBUG_END("NewtonEulerDS::computeJacobianqForces(double time) \n");
}

void NewtonEulerDS::computeJacobianvForces(double time)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianvForces(double time) \n");
  if (_jacobianvForces)
  {
    _jacobianvForces->zero();
    if (_jacobianFIntv)
    {
      computeJacobianFIntv(time);
      _jacobianvForces->setBlock(0,0,-1.0 * *_jacobianFIntv);
    }
    if (_jacobianMIntv)
    {
      computeJacobianMIntv(time);
      _jacobianvForces->setBlock(3,0,-1.0 * *_jacobianMIntv);
    }
    if (!nullifyFGyr)
    {
      if (_jacobianFGyrv)
      {
        //computeJacobianFGyrvByFD(time,_q,_v);
        computeJacobianFGyrv(time);
        *_jacobianvForces -= *_jacobianFGyrv;
      }
    }
    // std::cout << "_jacobianvForces : "<< std::endl;
    // _jacobianvForces->display();
  }
  //else nothing.
  DEBUG_END("NewtonEulerDS::computeJacobianvForces(double time) \n");
}

void NewtonEulerDS::computeJacobianFGyrv(double time)
{
  DEBUG_BEGIN("NewtonEulerDS::computeJacobianFGyrv(double time) \n");
  if (_jacobianFGyrv)
  {
    //Omega /\ I \Omega:
    _jacobianFGyrv->zero();
    SiconosVector omega(3);
    omega.setValue(0, _v->getValue(3));
    omega.setValue(1, _v->getValue(4));
    omega.setValue(2, _v->getValue(5));
    SiconosVector Iomega(3);
    prod(*_I, omega, Iomega, true);
    SiconosVector ei(3);
    SiconosVector Iei(3);
    SiconosVector ei_Iomega(3);
    SiconosVector omega_Iei(3);

    /*See equation of DevNotes.pdf, equation with label eq:NE_nablaFL1*/
    for (int i = 0; i < 3; i++)
    {
      ei.zero();
      ei.setValue(i, 1.0);
      prod(*_I, ei, Iei, true);
      cross_product(omega, Iei, omega_Iei);
      cross_product(ei, Iomega, ei_Iomega);
      for (int j = 0; j < 3; j++)
        _jacobianFGyrv->setValue(3 + j, 3 + i, ei_Iomega.getValue(j) + omega_Iei.getValue(j));
    }
    // Check if Jacobian is valid. Warning to the transpose operation in
    // _jacobianFGyrv->setValue(3 + j, 3 + i, ei_Iomega.getValue(j) + omega_Iei.getValue(j));
  }
  //else nothing.
  DEBUG_EXPR(_jacobianFGyrv->display());
  // _jacobianFGyrv->display();
  // SP::SimpleMatrix jacobianFGyrtmp (new SimpleMatrix(*_jacobianFGyrv));
  // computeJacobianFGyrvByFD(time, _q, _v);
  // jacobianFGyrtmp->display();
  // std::cout << "#################  " << (*jacobianFGyrtmp - *_jacobianFGyrv).normInf() << std::endl;
  // assert((*jacobianFGyrtmp - *_jacobianFGyrv).normInf()< 1e-10);
  DEBUG_END("NewtonEulerDS::computeJacobianFGyrv(double time) \n");
}


void NewtonEulerDS::display() const
{
  std::cout << "=====> NewtonEuler System display (number: " << _number << ")." <<std::endl;
  std::cout << "- _n : " << _n <<std::endl;
  std::cout << "- q " <<std::endl;
  if (_q) _q->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- q0 " <<std::endl;
  if (_q0) _q0->display();
  std::cout << "- v " <<std::endl;
  if (_v) _v->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- v0 " <<std::endl;
  if (_v0) _v0->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- dotq " <<std::endl;
  if (_dotq) _dotq->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[0] " <<std::endl;
  if (_p[0]) _p[0]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[1] " <<std::endl;
  if (_p[1]) _p[1]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "- p[2] " <<std::endl;
  if (_p[2]) _p[2]->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "mass :" <<  _scalarMass <<std::endl;
  std::cout << "Inertia :" <<std::endl;
  if (_I) _I->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "===================================== " <<std::endl;
}

// --- Functions for memory handling ---
void NewtonEulerDS::initMemory(unsigned int steps)
{
  DynamicalSystem::initMemory(steps);

  if (steps == 0)
    std::cout << "Warning : NewtonEulerDS::initMemory with size equal to zero" <<std::endl;
  else
  {
    _qMemory.reset(new SiconosMemory(steps, _qDim));
    _vMemory.reset(new SiconosMemory(steps, _n));
    _forcesMemory.reset(new SiconosMemory(steps, _n));
    _dotqMemory.reset(new SiconosMemory(steps, _qDim));
    swapInMemory();
  }
}

void NewtonEulerDS::swapInMemory()
{
  //  _xMemory->swap(_x[0]);
  _qMemory->swap(*_q);
  _vMemory->swap(*_v);
  _dotqMemory->swap(*_dotq);
  _forcesMemory->swap(*_forces);
}

void NewtonEulerDS::resetAllNonSmoothPart()
{
  if (_p[1])
    _p[1]->zero();
  else
    _p[1].reset(new SiconosVector(_n));
}
void NewtonEulerDS::resetNonSmoothPart(unsigned int level)
{
  if (_p[level]->size() > 0)
    _p[level]->zero();
}


void NewtonEulerDS::computeT()
{
  ::computeT(_q,_T);
}



void NewtonEulerDS::computeTdot()
{
  if (!_Tdot)
  {
    _Tdot.reset(new SimpleMatrix(_qDim, _n));
    _Tdot->zero();
  }

  ::computeT(_dotq,_Tdot);
}


void NewtonEulerDS::normalizeq()
{
  ::normalizeq(_q);
}




void NewtonEulerDS::setComputeJacobianFIntqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntqPtr, pluginPath,functionName);
  _pluginJacqFInt->setComputeFunction(pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianFIntvFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntvPtr, pluginPath,functionName);
  _pluginJacvFInt->setComputeFunction(pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianFIntqFunction(FInt_NE fct)
{
  _pluginJacqFInt->setComputeFunction((void *)fct);
}
void NewtonEulerDS::setComputeJacobianFIntvFunction(FInt_NE fct)
{
  _pluginJacvFInt->setComputeFunction((void *)fct);
}


void NewtonEulerDS::setComputeJacobianMIntqFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntqPtr, pluginPath,functionName);
  _pluginJacqMInt->setComputeFunction(pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianMIntvFunction(const std::string&  pluginPath, const std::string&  functionName)
{
  //    Plugin::setFunction(&computeJacobianFIntvPtr, pluginPath,functionName);
  _pluginJacvMInt->setComputeFunction(pluginPath, functionName);
}
void NewtonEulerDS::setComputeJacobianMIntqFunction(FInt_NE fct)
{
  _pluginJacqMInt->setComputeFunction((void *)fct);
}
void NewtonEulerDS::setComputeJacobianMIntvFunction(FInt_NE fct)
{
  _pluginJacvMInt->setComputeFunction((void *)fct);
}


double NewtonEulerDS::computeKineticEnergy()
{
  DEBUG_BEGIN("NewtonEulerDS::computeKineticEnergy()\n");
  assert(_v);
  assert(_massMatrix);
  DEBUG_EXPR(_v->display());
  DEBUG_EXPR(_massMatrix->display());

  SiconosVector tmp(6);
  prod(*_massMatrix, *_v, tmp, true);
  double K =0.5*inner_prod(tmp,*_v);

  DEBUG_PRINTF("Kinetic Energy = %e\n", K);
  DEBUG_END("NewtonEulerDS::computeKineticEnergy()\n");
  return K;
}
void NewtonEulerDS::setBoundaryConditions(SP::BoundaryCondition newbd)
{
  if(!_boundaryConditions){
    std::cout << "Warning : NewtonEulerDS::setBoundaryConditions. old boundary conditions were pre-existing" <<std::endl;
  }
  _boundaryConditions = newbd;
};
