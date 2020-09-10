#include "OccR.hpp"
#include "ContactPoint.hpp"
#include "OccContactShape.hpp"
#include "ContactShapeDistance.hpp"
#include "WhichGeometer.hpp"
#include "RuntimeException.hpp"
#include <limits>
#include <iostream>
#include <boost/typeof/typeof.hpp>

// #define  DEBUG_MESSAGES
#include "debug.h"

OccR::OccR(const ContactPoint& contact1,
           const ContactPoint& contact2,
           const DistanceCalculatorType& distance_calculator) :
  NewtonEuler3DR(),
  _contact1(contact1),
  _contact2(contact2),
  _geometer(),
  _offset1(0.),
  _offset2(0.)
{
  DEBUG_BEGIN("OccR::OccR(const ContactPoint& contact1, const ContactPoint& contact2,                         const DistanceCalculatorType& distance_calculator)\n");
  switch(Type::value(distance_calculator))
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
    
  DEBUG_END("OccR::OccR(const ContactPoint& contact1, const ContactPoint& contact2,                         const DistanceCalculatorType& distance_calculator)\n");
}


void OccR::computeh(double time, const BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("OccR::computeh(double time, BlockVector& q0, SiconosVector& y)\n");
  this->_contact2.contactShape().accept(*this->_geometer);

  
  
  ContactShapeDistance& dist = this->_geometer->answer;

  DEBUG_PRINTF("---->%g P1=(%g, %g, %g) P2=(%g,%g,%g) N=(%g, %g, %g)\n", dist.value,
               dist.x1, dist.y1, dist.z1,
               dist.x2, dist.y2, dist.z2,
               dist.nx, dist.ny, dist.nz);

  _Pc1->setValue(0, dist.x1 + _offset1*dist.nx);
  _Pc1->setValue(1, dist.y1 + _offset1*dist.ny);
  _Pc1->setValue(2, dist.z1 + _offset1*dist.nz);
  _Pc2->setValue(0, dist.x2 - _offset2*dist.nx);
  _Pc2->setValue(1, dist.y2 - _offset2*dist.ny);
  _Pc2->setValue(2, dist.z2 - _offset2*dist.nz);

  _Nc->setValue(0, dist.nx);
  _Nc->setValue(1, dist.ny);
  _Nc->setValue(2, dist.nz);

  dist.value -= (_offset1+_offset2);

  y.setValue(0, dist.value);

  DEBUG_EXPR(y.display(););
  DEBUG_EXPR(_Nc->display(););
  DEBUG_END("OccR::computeh(double time, BlockVector& q0, SiconosVector& y)\n");

}
