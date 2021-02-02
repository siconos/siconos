// This is the implementation for the class <MultipleImpactNSL>
#include "MultipleImpactNSL.hpp"
#include <iostream>

MultipleImpactNSL::MultipleImpactNSL(): NonSmoothLaw(1) {}
// Constructor with parameters
MultipleImpactNSL::MultipleImpactNSL(double newResCof, double newStiff, double newElasCoeff, unsigned int _dim): NonSmoothLaw(_dim)
{
  _ResCof = newResCof;
  _Stiff = newStiff;
  _ElasCof = newElasCoeff;
  // Throw the exceptions
  if((_ResCof < 0.0) || (_ResCof > 1.0))
    THROW_EXCEPTION("In MultipleImpactNSL, the restitution coefficient must be between 0.0 and 1.0!");
  if(_Stiff < 0.0)
    THROW_EXCEPTION("In MultipleImpactNSL, the stiffness must be positive!");
  if(_ElasCof < 0.0)
    THROW_EXCEPTION("In MultipleImpactNSL, the elasticity coefficient must be positive!");
}
// Destructor
MultipleImpactNSL::~MultipleImpactNSL() {}
//
void MultipleImpactNSL::setResCof(double newResCof)
{
  _ResCof = newResCof;
  if((_ResCof < 0.0) || (_ResCof > 1.0))
    THROW_EXCEPTION("MultipleImpactNSL::setResCof, the restitution coefficient must be between 0.0 and 1.0!");
}
//
void MultipleImpactNSL::setStiff(double newStiff)
{
  _Stiff = newStiff;
  if(_Stiff < 0.0)
    THROW_EXCEPTION("MultipleImpactNSL::setStiff, the stiffness must be positive!");
}
//
void MultipleImpactNSL::setElasCoeff(double _newElasCoef)
{
  _ElasCof = _newElasCoef;
  if(_newElasCoef < 0.0)
    THROW_EXCEPTION("MultipleImpactNSL::setElasCoeff, the elasticity coefficient must be positive!");
}
//
void MultipleImpactNSL::display() const
{
  std::cout << "===============================MultipleImpactNSL====================================" <<std::endl;
  std::cout << "Value of the energytical restitution coefficient at contact is :" << _ResCof <<std::endl;
  std::cout << "Value of the stiffness at contact is :"  << _Stiff <<std::endl;
  std::cout << "====================================================================================" <<std::endl;
}
bool MultipleImpactNSL::isVerified() const
{
  bool res = false;
  THROW_EXCEPTION("MultipleImpactNSL::isVerified is not yet implemented!");
  return res;
}
