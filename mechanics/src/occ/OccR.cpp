#include "OccR.hpp"
#include "ContactPoint.hpp"
#include "OccContactShape.hpp"
#include "ContactShapeDistance.hpp"
#include "WhichGeometer.hpp"
#include "RuntimeException.hpp"
#include <limits>
#include <iostream>
#include <boost/typeof/typeof.hpp>

OccR::OccR(const ContactPoint& contact1,
           const ContactPoint& contact2,
           const DistanceCalculatorType& distance_calculator) :
  NewtonEulerFrom3DLocalFrameR(),
  _contact1(contact1),
  _contact2(contact2),
  _geometer(),
  _normalFromFace1(false),
  _offsetp1(false),
  _offset(0.1)
{
  switch (Type::value(distance_calculator))
  {
  case Type::OccDistanceType:
    this->_geometer = ask<WhichGeometer<OccDistanceType> >(contact1.contactShape());
    break;
  case Type::CadmbtbDistanceType:
    this->_geometer = ask<WhichGeometer<CadmbtbDistanceType> >(contact1.contactShape());
    break;
  default:
    RuntimeException::selfThrow("OccR: Unknown distance calculator");
  }
  this->_contact2.contactShape().accept(*this->_geometer);
}


void OccR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  this->_contact2.contactShape().accept(*this->_geometer);

  ContactShapeDistance& dist = this->_geometer->answer;

  printf("---->%g P1=(%g, %g, %g) P2=(%g,%g,%g) N=(%g, %g, %g)\n", dist.value,
          dist.x1, dist.y1, dist.z1,
          dist.x2, dist.y2, dist.z2,
          dist.nx, dist.ny, dist.nz);

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
  _Nc->setValue(0, n1x);
  _Nc->setValue(1, n1y);
  _Nc->setValue(2, n1z);

  dist.value -= _offset;

  y.setValue(0, dist.value);

}
