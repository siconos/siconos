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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 *
 */

// for M_PI
#define _USE_MATH_DEFINES
#include <math.h>

#include "SphereLDS.hpp"

#define _2PI  2.0*M_PI

void normalize(SP::SiconosVector q, unsigned int i);

void normalize(SP::SiconosVector q, unsigned int i)
{

  q->setValue(i, fmod(q->getValue(i), _2PI));

  assert(fabs(q->getValue(i)) - std::numeric_limits<double>::epsilon() >= 0.);
  assert(fabs(q->getValue(i)) < _2PI);

}

void SphereLDS::computeMass()
{

  SP::SiconosVector qold;


  if (qMemory() && qMemory()->nbVectorsInMemory() >= 1)
    qold = qMemory()->getSiconosVector(0);

  normalize(q(), 3);
  normalize(q(), 4);
  normalize(q(), 5);

  if (qold)
  {
    normalize(qold, 3);
    normalize(qold, 4);
    normalize(qold, 5);
  }

  double theta = q()->getValue(3);

  assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);
  //assert (fabs(theta) - _2PI < 0.);

  (*_mass)(4, 5) = (*_mass)(5, 4) = I * cos(theta);

}

void SphereLDS::computeFGyr()
{
  SphereLDS::computeFGyr(q(), velocity());
}

void SphereLDS::computeFGyr(SP::SiconosVector q, SP::SiconosVector v)
{

  assert(q->size() == 6);
  assert(v->size() == 6);

  //  normalize(q,3);
  //normalize(q,4);
  //normalize(q,5);

  double theta    = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot   = v->getValue(4);
  double psidot   = v->getValue(5);

  double sintheta   = sin(theta);

  (*_fGyr)(0) = (*_fGyr)(1) = (*_fGyr)(2) = 0;

  (*_fGyr)(3) = I * psidot * phidot * sintheta;
  (*_fGyr)(4) = -I * psidot * thetadot * sintheta;
  (*_fGyr)(5) = -I * phidot * thetadot * sintheta;
}



void SphereLDS::computeJacobianFGyrq()
{

  SphereLDS::computeJacobianFGyrq(q(), velocity());
}
void SphereLDS::computeJacobianFGyrqDot()
{

  SphereLDS::computeJacobianFGyrqDot(q(), velocity());
}

void SphereLDS::computeJacobianFGyrq(SP::SiconosVector q, SP::SiconosVector v)
{
  double theta    = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot   = v->getValue(4);
  double psidot   = v->getValue(5);

  double costheta = cos(theta);

  _jacobianFGyrq->zero();

  (*_jacobianFGyrq)(3, 3) = -I * psidot * phidot * costheta;
  (*_jacobianFGyrq)(4, 3) = I * psidot * thetadot * costheta;
  (*_jacobianFGyrq)(5, 3) = I * psidot * thetadot * costheta;


}
void SphereLDS::computeJacobianFGyrqDot(SP::SiconosVector q, SP::SiconosVector v)
{
  double theta    = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot   = v->getValue(4);
  double psidot   = v->getValue(5);

  double sintheta   = sin(theta);

  _jacobianFGyrqDot->zero();


  (*_jacobianFGyrqDot)(3, 3) = 0;
  (*_jacobianFGyrqDot)(3, 4) = I * psidot * sintheta;
  (*_jacobianFGyrqDot)(3, 5) = I * phidot * sintheta;

  (*_jacobianFGyrqDot)(4, 3) = -I * psidot * sintheta;
  (*_jacobianFGyrqDot)(4, 4) = 0;
  (*_jacobianFGyrqDot)(4, 5) = -I * thetadot * sintheta;

  (*_jacobianFGyrqDot)(5, 3) =  -I * phidot * sintheta;
  (*_jacobianFGyrqDot)(5, 4) =  -I * thetadot * sintheta;
  (*_jacobianFGyrqDot)(5, 5) = 0;

}


SphereLDS::SphereLDS(double r, double m,
                     SP::SiconosVector qinit,
                     SP::SiconosVector vinit)
  : LagrangianDS(qinit, vinit), radius(r), massValue(m)
{

  normalize(q(), 3);
  normalize(q(), 4);
  normalize(q(), 5);
  _ndof = 6;

  assert(qinit->size() == _ndof);
  assert(vinit->size() == _ndof);

  _mass.reset(new SimpleMatrix(_ndof, _ndof));
  _mass->zero();
  I = massValue * radius * radius * 2. / 5.;
  (*_mass)(0, 0) = (*_mass)(1, 1) = (*_mass)(2, 2) = massValue;    ;
  (*_mass)(3, 3) = (*_mass)(4, 4) = (*_mass)(5, 5) = I;

  computeMass();

  _jacobianFGyrq.reset(new SimpleMatrix(_ndof, _ndof));
  _jacobianFGyrqDot.reset(new SimpleMatrix(_ndof, _ndof));

  _fGyr.reset(new SiconosVector(_ndof));
  _fGyr->zero();

  computeFGyr();


}

SphereLDS::~SphereLDS()
{}
