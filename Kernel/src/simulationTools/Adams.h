#ifndef ADAMS_H
#define ADAMS_H

#include "OneStepIntegrator.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "AdamsXML.h"

//using namespace std;

/** \class Adams
 *  \brief Adams is a kind of multi-step integrator.
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class Adams : public OneStepIntegrator
{
public:
  /** \fn Adams()
   *  \brief default constructor
   */
  Adams();

  /** \fn Adams(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor with XML object of the Adams
   *  \param OneStepIntegratorXML* : the XML object corresponding
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Adams(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem*);

  ~Adams();

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
  inline void setR(int r)
  {
    this->r = r;
  }


  /** \fn void saveIntegratorToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   */
  void saveIntegratorToXML();

  /** \fn void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*, DynamicalSystem*)
   *  \brief allows to create the Integrator Adams with an xml file, or the needed data
   *  \param OneStepNSProblemXML * : the XML object for this OneStepIntegrator
   *  \param TimeDiscretisation * : The NSDS which contains this OneStepIntegrator
   *  \exception RuntimeException
   */
  void createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation* td, DynamicalSystem* ds);//, Strategy * strategy = NULL);

  /** \fn void initialize()
  *  \brief initialization of the Adams integrator
  */
  void initialize();

  /** \fn void fillIntegratorWithIntegratorXML()
   *  \brief uses the OneStepIntegratorXML of the Adams Integrator to fill the fields of this Integrator
   *  \exception RuntimeException
   */
  void fillIntegratorWithIntegratorXML();

  /** \fn Adams* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Adams* convert(OneStepIntegrator* osi);

private:
  /**  */
  int r;

};

#endif // ADAMS_H
