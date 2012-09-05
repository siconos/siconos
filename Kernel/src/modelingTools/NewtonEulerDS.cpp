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
#include "NewtonEulerDS.hpp"
#include "BlockVector.hpp"
#include "BlockMatrix.hpp"
#include <boost/math/quaternion.hpp>
using namespace std;

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
NewtonEulerDS::NewtonEulerDS(): DynamicalSystem(6)
{
  _p.resize(3);
  zeroPlugin();
  //assert(0);
  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --


  _qDim = 7;
  _n = 6;


  // Current state
  _q.reset(new SiconosVector(_qDim));
  _deltaq.reset(new SiconosVector(_qDim));
  _v.reset(new SiconosVector(_n));

  _dotq.reset(new SiconosVector(_qDim));
  _residuFree.reset(new SiconosVector(_n));
  _massMatrix.reset(new SimpleMatrix(_n, _n));
  _luW.reset(new SimpleMatrix(_n, _n));
  _massMatrix->zero();
  _T.reset(new SimpleMatrix(_qDim, _n));

}
void NewtonEulerDS::internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix)
{
  _p.resize(3);
  zeroPlugin();
  // --- NEWTONEULER INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --

  _mass = mass;
  _qDim = 7;
  _n = 6;

  // Initial conditions
  _q0 = Q0;
  _v0 = Velocity0;

  _MObjToAbs.reset(new SimpleMatrix(3, 3));

  // Current state
  _q.reset(new SiconosVector(_qDim));
  _deltaq.reset(new SiconosVector(_qDim));
  _v.reset(new SiconosVector(_n));
  (*_q) = (*_q0);
  _dotq.reset(new SiconosVector(_qDim));
  _residuFree.reset(new SiconosVector(_n));
  _massMatrix.reset(new SimpleMatrix(_n, _n));
  _jacobianvFL.reset(new SimpleMatrix(_n, _n));
  _luW.reset(new SimpleMatrix(_n, _n));
  _massMatrix->zero();
  _massMatrix->setValue(0, 0, _mass);
  _massMatrix->setValue(1, 1, _mass);
  _massMatrix->setValue(2, 2, _mass);
  _I = inertialMatrix;
  Index dimIndex(2);
  dimIndex[0] = 3;
  dimIndex[1] = 3;
  Index startIndex(4);
  startIndex[0] = 0;
  startIndex[1] = 0;
  startIndex[2] = 3;
  startIndex[3] = 3;
  setBlock(_I, _massMatrix, dimIndex, startIndex);

  _T.reset(new SimpleMatrix(_qDim, _n));
  _T->zero();
  _T->setValue(0, 0, 1.0);
  _T->setValue(1, 1, 1.0);
  _T->setValue(2, 2, 1.0);
  updateT();
  initForces();
}
NewtonEulerDS::NewtonEulerDS(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double  mass, SP::SiconosMatrix inertialMatrix):
  DynamicalSystem(6)
{
  internalInit(Q0, Velocity0, mass, inertialMatrix);
}

void NewtonEulerDS::zeroPlugin()
{
  computeJacobianFIntqPtr = NULL;
  computeJacobianFIntqDotPtr = NULL;
  _pluginFExt.reset(new PluggedObject());
  _pluginMExt.reset(new PluggedObject());
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

  // q0 and velocity0
  if (! _q0 || ! _v0)
  {
    RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - initial conditions are badly set.");
    output = false;
  }


  // fInt
  //   if( ( _fInt && computeFIntPtr) && ( ! _jacobianFIntq || ! _jacobianFIntqDot ) )
  //     // ie if fInt is defined and not constant => its Jacobian must be defined (but not necessarily plugged)
  //     {
  //       RuntimeException::selfThrow("NewtonEulerDS::checkDynamicalSystem - You defined fInt but not its Jacobian (according to q and velocity).");
  //       output = false;
  //     }


  if (!output) cout << "NewtonEulerDS Warning: your dynamical system seems to be uncomplete (check = false)" << endl;
  return output;
}

void NewtonEulerDS::initializeNonSmoothInput(unsigned int level)
{

  if (level == 0)
  {
    if (!_p[0])
      _p[0].reset(new SiconosVector(_qDim));
  }
  else
  {
    if (!_p[level])
      _p[level].reset(new SiconosVector(_n));
  }
}

void NewtonEulerDS::initForces()
{

  _forces.reset(new SiconosVector(_n));

  _jacobianvFL.reset(new SimpleMatrix(_n, _qDim));
  _jacobianqDotForces.reset(new SimpleMatrix(_n, _qDim));
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

  //   if ( computeFIntPtr && ! _fInt)
  //     _fInt.reset(new SiconosVector(_n));
  //   if (computeJacobianFIntqPtr && !_jacobianFIntq)
  //     _jacobianFIntq.reset(new SimpleMatrix(_n,_n));
  //   if (computeJacobianFIntqDotPtr && !_jacobianFIntqDot)
  //     _jacobianFIntqDot.reset(new SimpleMatrix(_n,_n));


  //
  if (!_workFree)
    _workFree.reset(new SiconosVector(getDim()));
  // Memory allocation for fL and its jacobians.

  // Set links to variables of top-class DynamicalSystem.
  // Warning: this results only in pointers links. No more memory allocation for vectors or matrices.
  connectToDS(); // note that connection can not be done during constructor call, since user can complete the ds after (add plugin or anything else).
  checkDynamicalSystem();

  // Initialize memory vectors
  initMemory(sizeOfMemory);

  initRhs(time);

}



void NewtonEulerDS::computeFExt(const double time)
{
  if (_pluginFExt->fPtr)
    ((Fext)_pluginFExt->fPtr)(time, &(*_q)(0), &(*_fExt)(0),  &(*_q0)(0));
}

void NewtonEulerDS::computeMExt(const double time)
{
  if (_pluginMExt->fPtr)
    ((Fext)_pluginMExt->fPtr)(time, &(*_q)(0), &(*_mExt)(0),  &(*_q0)(0));
}

void NewtonEulerDS::computeRhs(const double time, const bool isDSup)
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
      computeMExt(time);
      SiconosVector aux(3);
      updateMObjToAbs();
      prod(*_mExt, *_MObjToAbs, aux);
      *_mExt = aux;
      _forces->setBlock(3, *_mExt);
    }
    /*computation of \Omega vectortiel I \Omega*/
    if (_I)
    {
      // printf("NewtonEulerDS::computeFL _I:\n");
      // _I->display();
      // _v->display();

      SiconosVector bufOmega(3);
      SiconosVector bufIOmega(3);
      SiconosVector buf(3);
      bufOmega.setValue(0, _v->getValue(3));
      bufOmega.setValue(1, _v->getValue(4));
      bufOmega.setValue(2, _v->getValue(5));
      prod(*_I, bufOmega, bufIOmega, true);
      cross_product(bufOmega, bufIOmega, buf);
      _forces->setValue(3, _forces->getValue(3) - buf.getValue(0));
      _forces->setValue(4, _forces->getValue(4) - buf.getValue(1));
      _forces->setValue(5, _forces->getValue(5) - buf.getValue(2));
    }
  }
  // else nothing.
}
void NewtonEulerDS::computeForces(double time, SP::SiconosVector q2, SP::SiconosVector v2)
{
  if (_forces)
  {
    _forces->zero();
    // 1 - Computes the required functions
    if (_fExt)
    {
      //      computeFExt(time);
      _forces->setBlock(0, *_fExt);
    }
    if (_mExt)
    {
      //      computeMExt(time);
      SiconosVector aux(3);
      updateMObjToAbs();
      prod(*_mExt, *_MObjToAbs, aux);
      *_mExt = aux;
      _forces->setBlock(3, *_mExt);
    }
    /*computation of \Omega vectortiel I \Omega*/
    if (_I)
    {
      SiconosVector bufOmega(3);
      SiconosVector bufIOmega(3);
      SiconosVector buf(3);
      bufOmega.setValue(0, v2->getValue(3));
      bufOmega.setValue(1, v2->getValue(4));
      bufOmega.setValue(2, v2->getValue(5));
      prod(*_I, bufOmega, bufIOmega, true);
      cross_product(bufOmega, bufIOmega, buf);
      _forces->setValue(3, _forces->getValue(3) - buf.getValue(0));
      _forces->setValue(4, _forces->getValue(4) - buf.getValue(1));
      _forces->setValue(5, _forces->getValue(5) - buf.getValue(2));
    }
  }
}


void NewtonEulerDS::computeJacobianvFL(double time)
{
  if (_jacobianvFL)
  {
    //Omega /\ I \Omega:
    _jacobianvFL->zero();
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
        _jacobianvFL->setValue(3 + i, 3 + j, ei_Iomega.getValue(j) + omega_Iei.getValue(j));
    }
  }
  //else nothing.
}
void NewtonEulerDS::computeJacobianqDotForces(double time)
{
  if (_jacobianqDotForces)
  {
    //      computeJacobianFIntqDot(time);

    // not true!
    // if( jacobianFL[i].use_count() == 1 )
    {
      //if not that means that jacobianFL[i] is already (pointer-)connected with
      // either jacobianFInt
      _jacobianqDotForces->zero();
      //    if( _jacobianFIntqDot )
      //      *_jacobianqDotForces-=*_jacobianFIntqDot;
    }
  }
  //else nothing.
}
// void NewtonEulerDS::computeJacobianZFL( double time){
//    RuntimeException::selfThrow("NewtonEulerDS::computeJacobianZFL - not implemented");
// }

void NewtonEulerDS::saveSpecificDataToXML()
{

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
    _forcesMemory.reset(new SiconosMemory(steps));
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
  _forcesMemory->swap(_forces);

}

NewtonEulerDS* NewtonEulerDS::convert(DynamicalSystem* ds)
{
  NewtonEulerDS* lnlds = dynamic_cast<NewtonEulerDS*>(ds);
  return lnlds;
}

void NewtonEulerDS::resetNonSmoothPart()
{
  if (_p[1])
    _p[1]->zero();
  else
    _p[1].reset(new SiconosVector(_n));
}
void NewtonEulerDS::resetNonSmoothPart(unsigned int level)
{
  if (_p[level])
    _p[level]->zero();
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
void NewtonEulerDS::normalizeq()
{
  double normq = sqrt(_q->getValue(3) * _q->getValue(3) + _q->getValue(4) * _q->getValue(4) + _q->getValue(5) * _q->getValue(5) + _q->getValue(6) * _q->getValue(6));
  assert(normq > 0);
  normq = 1 / normq;
  _q->setValue(3, _q->getValue(3) * normq);
  _q->setValue(4, _q->getValue(4) * normq);
  _q->setValue(5, _q->getValue(5) * normq);
  _q->setValue(6, _q->getValue(6) * normq);
}
void NewtonEulerDS::updateMObjToAbs()
{
  double q0 = _q->getValue(3);
  double q1 = _q->getValue(4);
  double q2 = _q->getValue(5);
  double q3 = _q->getValue(6);

  ::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
  ::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
  ::boost::math::quaternion<double>    quatx(0, 1, 0, 0);
  ::boost::math::quaternion<double>    quaty(0, 0, 1, 0);
  ::boost::math::quaternion<double>    quatz(0, 0, 0, 1);
  ::boost::math::quaternion<double>    quatBuff;
  /*See equation with label eq:newton_Mobjtoabs from the DevNote.pdf, chapter Gradiant computaion, case oif NewtonEuler with quaternion*/
  quatBuff = quatcQ * quatx * quatQ;
  _MObjToAbs->setValue(0, 0, quatBuff.R_component_2());
  _MObjToAbs->setValue(0, 1, quatBuff.R_component_3());
  _MObjToAbs->setValue(0, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quaty * quatQ;
  _MObjToAbs->setValue(1, 0, quatBuff.R_component_2());
  _MObjToAbs->setValue(1, 1, quatBuff.R_component_3());
  _MObjToAbs->setValue(1, 2, quatBuff.R_component_4());
  quatBuff = quatcQ * quatz * quatQ;
  _MObjToAbs->setValue(2, 0, quatBuff.R_component_2());
  _MObjToAbs->setValue(2, 1, quatBuff.R_component_3());
  _MObjToAbs->setValue(2, 2, quatBuff.R_component_4());

}
