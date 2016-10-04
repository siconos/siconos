#include "OccR.hpp"
#include "ContactPoint.hpp"
#include "OccContactShape.hpp"
#include "ContactShapeDistance.hpp"
#include "WhichGeometer.hpp"
#include <limits>
#include <iostream>
#include <boost/typeof/typeof.hpp>

OccR::OccR(const ContactPoint& contact1,
           const ContactPoint& contact2,
           const DistanceCalculatorType& distance_calculator) :
  NewtonEulerFrom3DLocalFrameR(),
  _contact1(contact1),
  _contact2(contact2),
  _geometer(ask<WhichGeometer<BOOST_TYPEOF(distance_calculator)> >(contact1.contactShape())),
  _normalFromFace1(false),
  _offsetp1(false),
  _offset(0.002)
{
}


void OccR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  this->_contact2.contactShape().accept(*this->_geometer);

  ContactShapeDistance& dist = this->_geometer->answer;

  double& X1 = dist.x1;
  double& Y1 = dist.y1;
  double& Z1 = dist.z1;

  double& X2 = dist.x2;
  double& Y2 = dist.y2;
  double& Z2 = dist.z2;

  double& n1x = dist.nx;
  double& n1y = dist.ny;
  double& n1z = dist.nz;

  if(_offsetp1)
  {
    _Pc1->setValue(0, X1+_offset*n1x);
    _Pc1->setValue(1, Y1+_offset*n1y);
    _Pc1->setValue(2, Z1+_offset*n1z);
    _Pc2->setValue(0, X2);
    _Pc2->setValue(1, Y2);
    _Pc2->setValue(2, Z2);
  }
  else
  {
    _Pc1->setValue(0, X1);
    _Pc1->setValue(1, Y1);
    _Pc1->setValue(2, Z1);
    _Pc2->setValue(0, X2-_offset*n1x);
    _Pc2->setValue(1, Y2-_offset*n1y);
    _Pc2->setValue(2, Z2-_offset*n1z);
  }

  /* cf comments from O. Bonnefon */
  _Nc->setValue(0, -n1x);
  _Nc->setValue(1, -n1y);
  _Nc->setValue(2, -n1z);

  dist.value -= _offset;

  y.setValue(0, dist.value);

}
