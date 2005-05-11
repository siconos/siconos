#ifndef MOREAU_H
#define MOREAU_H

#include "OneStepIntegrator.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "MoreauXML.h"
#include "check.h"


const int MOREAUSTEPSINMEMORY = 1;

//using namespace std;

/** \class Moreau
 *  \brief It's a kind of single-step Integrator
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 * \todo Add the LU Factorization of W in the initialization,  and adapt the resolution in the iteration
 */
class Moreau : public OneStepIntegrator
{
public:

  /** \fn Moreau(OneStepIntegratorXML*,TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  /** \fn Moreau(TimeDiscretisation*, DynamicalSystem* , const double& theta)
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   *  \param Theta value
   */
  Moreau(TimeDiscretisation*, DynamicalSystem*, const double& theta);

  /** \fn ~Moreau()
   *  \brief destructor
   */
  ~Moreau();

  // --- GETTERS/SETTERS ---
  // -- W --

  /** \fn  const SiconosMatrix getW(void) const
   *  \brief get the value of W
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getW(void) const
  {
    return *(this->W);
  }

  /** \fn SiconosMatrix* getWPtr(void) const
   *  \brief get W
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getWPtr(void) const
  {
    return this->W;
  }

  /** \fn void setW (const SiconosMatrix& newValue)
   *  \brief set the value of W to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setW(const SiconosMatrix& newValue)
  {
    *(this->W) = newValue;
  }

  /** \fn void setWPtr(SiconosMatrix* newPtr)
   *  \brief set W to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setWPtr(SiconosMatrix *newPtr)
  {
    delete W;
    W = 0;
    W = newPtr;
  }

  // -- theta --

  /** \fn const double getTheta() const
   *  \brief allows to get the double value theta of the Moreau's Integrator
   *  \return double value theta
   */
  inline const double getTheta() const
  {
    return theta;
  }

  /** \fn double setTheta(const double&)
   *  \brief set the value of theta
   *  \param ref on a double
   */
  inline void setTheta(const double& newTheta)
  {
    theta = newTheta;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void initialize()
   *  \brief initialization of the Moreau integrator; for linear time invariant systems, we compute time invariant operator (example : W)
   *  \todo LU factorization of time invariant operator (example : W)
   */
  void initialize();

  /** \fn void computeFreeState()
   *  \brief integrates the Dynamical System linked to this integrator without boring the constraints
   */
  void computeFreeState();

  /** \fn void integrate()
   *  \brief makes computations to integrate the data of a Dynamical System with the Moreau Integrator
   */
  void integrate();

  /** \fn void updateState()
   *  \brief updates the state of the Dynamical System
   */
  void updateState();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveIntegratorToXML();

  /** \fn void saveWToXML()
   *  \brief copy the matrix W of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveWToXML();


  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the Moreau Integrator to fill the fields of this Integrator
   *  \exception RuntimeException
   */
  void fillIntegratorWithIntegratorXML();

  /** \fn display()
   *  \brief Displays the data of the Moreau's integrator
   */
  void display() const;

  /** \fn Moreau* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, 0 otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

private:
  /** \fn Moreau()
   *  \brief Default constructor
   */
  Moreau();

  /** a specific matrix of the Moreau Integrator */
  SiconosMatrix *W;

  /** parameter of the theta method */
  double theta;
};

#endif // MOREAU_H
