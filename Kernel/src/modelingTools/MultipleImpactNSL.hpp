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
  //Elasticity coefficient
  double _ElasCof;
public:
  // Default Constructor
  MultipleImpactNSL();
  // Constructor with parameters
  MultipleImpactNSL(double, double, double, unsigned int _dim = 1);
  // Destructor
  ~MultipleImpactNSL();
  // Get the value of the energytical restitution coefficientx
  inline double ResCof() const
  {
    return _ResCof;
  };
  // Get the value of the stiffness
  inline double Stiff() const
  {
    return _Stiff;
  };
  // Get the value of the elasticity coefficient
  inline double ElasCof() const
  {
    return _ElasCof;
  }
  // Set the value to the restitution coefficient
  void setResCof(double newResCof);
  // Set the value to the stiffness
  void setStiff(double newStiff);
  // Set the value to the elasticity cofficient
  void setElasCoeff(double _newElasCoef);
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














