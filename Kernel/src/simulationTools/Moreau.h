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

using namespace std;

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

  /** \fn Moreau()
   *  \brief Default constructor
   */
  Moreau();

  /** \fn Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the Moreau
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Moreau(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  ~Moreau();

  //getter/setter
  /** \fn SiconosMatrix getW(void)
   *  \brief allows to get the SiconosMatrix W of the Moreau's Integrator
   *  \return the SiconosMatrix W
   */
  inline SiconosMatrix getW(void) const
  {
    return this->W;
  };

  /** \fn SiconosMatrix* getWPtr(void)
   *  \brief allows to get the SiconosMatrix* W of the Moreau's Integrator
   *  \return the SiconosMatrix* W
   */
  SiconosMatrix* getWPtr(void);

  /** \fn void setW(SiconosMatrix)
   *  \brief allows to set the SiconosMatrix W of the Moreau's Integrator
   *  \param the SiconosMatrix to set W
   */
  inline void setW(const SiconosMatrix& W)
  {
    this->W = W;
  };

  /** \fn double getTheta(void)
   *  \brief allows to get the double value theta of the Moreau's Integrator
   *  \return double value theta
   */
  inline double getTheta(void) const
  {
    return this->theta;
  }

  /** \fn double getTheta(void)
   *  \brief allows to set the double value theta of the Moreau's Integrator
   *  \param double value to set
   */
  inline void setTheta(const double theta)
  {
    this->theta = theta;
  }


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

  /** \fn void createOneStepIntegrator(OneStepIntegratorXML * osiXML,
      TimeDiscretisation* td, DynamicalSystem* ds, double theta)
   *  \brief allows to create the Integrator Moreau with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepIntegrator
   *  \param double :  the value for theta of this OneStepIntegrator
   *  \param TimeDiscretisation * : The NSDS which contains this OneStepIntegrator
   *  \exception RuntimeException
   */
  void createOneStepIntegrator(OneStepIntegratorXML * osiXML,
                               TimeDiscretisation* td, DynamicalSystem* ds,
                               /*int r=-1,*/ double theta = -1.0); //, Strategy * strategy = NULL);

  /** \fn Moreau* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Moreau* convert(OneStepIntegrator* osi);

private:
  /** a specific matrix of the Moreau Integrator */
  SiconosMatrix W;

  /** parameter of the theta method */
  double theta;
};

#endif // MOREAU_H
