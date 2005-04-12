#ifndef Lsodar_H
#define Lsodar_H

#include "OneStepIntegrator.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "LsodarXML.h"
#include "check.h"

//#include "C_file_for_Lsodar.h"

//using namespace std;

/** \class Lsodar
 *  \brief It's a kind of single-step Integrator
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 */
class Lsodar : public OneStepIntegrator
{
public:

  /** \fn Lsodar()
   *  \brief default constructor
   */
  Lsodar();

  /** \fn Lsodar(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the Lsodar
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   Integrator
   */
  Lsodar(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  ~Lsodar();


  /** \fn double getR()
  *   \brief Return the r of the OneStepIntegrator
  *   \return int : the value of r
  */
  inline int getR(void) const
  {
    return this->r;
  }

  /** \fn void setR(int r)
  *   \brief Return the r of OneStepIntegrator
  *   \param double : the value to set r
  */
  inline void setR(const int r)
  {
    this->r = r;
  }

  /** \fn void computeFreeState()
  *   \brief compute the free state of the dynamical system
  */
  void computeFreeState();

  /** \fn void integrate()
  *   \brief integrates the dynamical system
  */
  void integrate();

  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveIntegratorToXML();

  /** \fn void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*, DynamicalSystem*)
   *  \brief allows to create the Integrator Lsodar with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepIntegrator
   *  \param TimeDiscretisation * : The NSDS which contains this OneStepIntegrator
   *  \exception RuntimeException
   */
  void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation* td, DynamicalSystem* ds);//, Strategy * strategy = NULL);

  /** \fn void initialize()
  *  \brief initialization of the Lsodar integrator
  */
  void initialize();

  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the Lsodar Integrator to fill the fields of this Integrator
   *  \exception RuntimeException
   */
  void fillIntegratorWithIntegratorXML();

  /** \fn Lsodar* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Lsodar* convert(OneStepIntegrator* osi);

};
typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);
extern "C" void tryfunction(fctPtr);

#endif // Lsodar_H
