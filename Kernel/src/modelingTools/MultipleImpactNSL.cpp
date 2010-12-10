// This is the implementation for the class <MultipleImpactNSL>
#include "MultipleImpactNSL.hpp"

MultipleImpactNSL::MultipleImpactNSL(): NonSmoothLaw(1), _ResCof(0.0), _Stiff(0.0) {}
// Constructor with parameters
MultipleImpactNSL::MultipleImpactNSL(double newResCof, double newStiff): NonSmoothLaw(1)
{
  _ResCof = newResCof;
  _Stiff = newStiff;
  // Throw the exceptions
  if ((_ResCof < 0.0) || (_ResCof > 1.0))
    RuntimeException::selfThrow("MultipleImpactNSL::_ResCof must be between 0.0 and 1.0!");
  if (_Stiff < 0.0)
    RuntimeException::selfThrow("MultipleImpactNSL::_Stiff must be positive!");
}
// Destructor
MultipleImpactNSL::~MultipleImpactNSL() {}
//
void MultipleImpactNSL::display() const
{
  cout << "===============================MultipleImpactNSL====================================" << endl;
  cout << "Value of the energytical restitution coefficient at contact is :" << _ResCof << endl;
  cout << "Value of the stiffness at contact is :"  << _Stiff << endl;
  cout << "====================================================================================" << endl;
}
bool MultipleImpactNSL::isVerified() const
{
  bool res = false;
  RuntimeException::selfThrow("MultipleImpactNSL::isVerified is not yet implemented!");
  return res;
}
//
void MultipleImpactNSL::saveNonSmoothLawToXML()
{
  RuntimeException::selfThrow("MultipleImpactNSL::saveNonSmoothLawToXML is not yet implemented!");
}
