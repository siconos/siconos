/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 *
 */

#include "Sphere.h"

static double _2PI = 2 * M_PI;


void normalize(SP::SiconosVector q, unsigned int i)
{

  q->setValue(i, fmod(q->getValue(i), _2PI));

  assert(fabs(q->getValue(i)) - std::numeric_limits<double>::epsilon() >= 0.);
  assert(fabs(q->getValue(i)) < _2PI);

}

void Sphere::computeMass()
{

  SP::SiconosVector q = getQPtr();
  SP::SiconosVector qold;


  if (getQMemoryPtr())
    qold = getQMemoryPtr()->getSiconosVector(0);

  normalize(q, 3);
  normalize(q, 4);
  normalize(q, 5);

  if (qold)
  {
    normalize(qold, 3);
    normalize(qold, 4);
    normalize(qold, 5);
  }

  double theta = getQPtr()->getValue(3);

  assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);
  //assert (fabs(theta) - _2PI < 0.);

  (*mass)(4, 5) = (*mass)(5, 4) = I * cos(theta);

}

void Sphere::computeNNL()
{
  Sphere::computeNNL(getQPtr(), getVelocityPtr());
}

void Sphere::computeNNL(SP::SiconosVector q, SP::SiconosVector v)
{

  assert(q->size() == 6);
  assert(v->size() == 6);

  //  normalize(q,3);
  //normalize(q,4);
  //normalize(q,5);

  double theta    = q->getValue(3);
  double phi      = q->getValue(4);
  double psi      = q->getValue(5);

  double thetadot = v->getValue(3);
  double phidot   = v->getValue(4);
  double psidot   = v->getValue(5);


  (*NNL)(0) = (*NNL)(1) = (*NNL)(2) = 0;

  (*NNL)(3) = I * psidot * phidot * sin(theta);
  (*NNL)(4) = -I * psidot * thetadot * sin(theta);
  (*NNL)(5) = -I * phidot * thetadot * sin(theta);
}

Sphere::Sphere(double r, double m,
               const SiconosVector& qinit,
               const SiconosVector& vinit)
  : LagrangianDS(qinit, vinit), radius(r), massValue(m)
{

  normalize(getQPtr(), 3);
  normalize(getQPtr(), 4);
  normalize(getQPtr(), 5);
  ndof = 6;

  assert(qinit.size() == ndof);
  assert(vinit.size() == ndof);

  double theta    = qinit.getValue(3);
  double phi      = qinit.getValue(4);
  double psi      = qinit.getValue(5);

  double thetadot = vinit.getValue(3);
  double phidot   = vinit.getValue(4);
  double psidot   = vinit.getValue(5);


  mass.reset(new PMMass(ndof, ndof));
  mass->zero();
  I = massValue * radius * radius * 2. / 5.;
  (*mass)(0, 0) = (*mass)(1, 1) = (*mass)(2, 2) = massValue;    ;
  (*mass)(3, 3) = (*mass)(4, 4) = (*mass)(5, 5) = I;

  computeMass();


  NNL.reset(new PVNNL(ndof));
  NNL->zero();

  computeNNL();


}

Sphere::~Sphere()
{}
