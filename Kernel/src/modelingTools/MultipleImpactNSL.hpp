// This is the header file to define the parameters of the multiple impact law
#ifndef _MULTIPLEIMPACTNSL_
#define _MULTIPLEIMPACTNSL_
//==================================================================================================================
#include "NonSmoothLaw.hpp"
#include <iostream>
using namespace std;
//==================================================================================================================
class MultipleImpactNSL : public NonSmoothLaw
{
private:
  //Energytical restitution coefficient
  double _ResCof;
  //Normal stiffness at contact
  double _Stiff;
public:
  // Default Constructor
  MultipleImpactNSL();
  // Constructor with parameters
  MultipleImpactNSL(double, double);
  // Destructor
  ~MultipleImpactNSL();
  // Get the value of the energytical restitution coefficient
  inline double ResCof() const
  {
    return _ResCof;
  };
  // Get the value of the stiffness
  inline double Stiff() const
  {
    return _Stiff;
  };
  // Set the value to the restitution coefficient
  inline void setResCof(double newResCof)
  {
    _ResCof = newResCof;
    if ((_ResCof < 0.0) || (_ResCof > 1.0))
      RuntimeException::selfThrow("MultipleImpactNSL::_ResCof must be between 0.0 and 1.0!");
  };
  // Set the value to the stiffness
  inline void setStiff(double newStiff)
  {
    _Stiff = newStiff;
    if (_Stiff < 0.0)
      RuntimeException::selfThrow("MultipleImpactNSL::_Stiff must be positive!");
  };
  //
  bool isVerified() const;
  //
  void saveNonSmoothLawToXML();
  // Display the information about the multiple impact law
  void display() const;
  //Visitors hook
  ACCEPT_STD_VISITORS();
};
DEFINE_SPTR(MultipleImpactNSL);
#endif














