#ifndef Lsodar_H
#define Lsodar_H

#include "OneStepIntegrator.h"
#include "LsodarXML.h"

/** \class Lsodar
 *  \brief It's a kind of single-step Integrator
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 */
class Lsodar : public OneStepIntegrator
{
public:


  /** \fn Lsodar(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object
   */
  Lsodar(OneStepIntegratorXML*);

  /** \fn Lsodar(TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   Integrator
   */
  Lsodar(TimeDiscretisation*, DynamicalSystem*);

  ~Lsodar();

  /** \fn void computeFreeState()
  *   \brief compute the free state of the dynamical system
  */
  void computeFreeState();

  /** \fn void integrate()
  *   \brief integrates the dynamical system
  */
  void integrate();

  /** \fn Lsodar* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Lsodar* convert(OneStepIntegrator* osi);

private:
  /** \fn Lsodar()
   *  \brief default constructor
   */
  Lsodar();
};


typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);
extern "C" void tryfunction(fctPtr);

#endif // Lsodar_H
