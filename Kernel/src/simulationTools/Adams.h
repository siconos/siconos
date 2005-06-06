#ifndef ADAMS_H
#define ADAMS_H

#include "OneStepIntegrator.h"
#include "AdamsXML.h"

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

  /** \fn Adams(OneStepIntegratorXML*)
   *  \brief constructor from xml file
   *  \param OneStepIntegratorXML* : the XML object corresponding
   */
  Adams(OneStepIntegratorXML*);
  /** \fn Adams(OneStepIntegratorXML*, TimeDiscretisation*, DynamicalSystem* )
   *  \brief constructor from a minimum set of data
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem linked to the OneStepIntegrator
   */
  Adams(TimeDiscretisation*, DynamicalSystem*);

  ~Adams();

  /** \fn double const getR() const
   *   \brief Return the r of the OneStepIntegrator
   *   \return int : the value of r
   */
  inline const int getR() const
  {
    return r;
  }

  /** \fn void setR(const int& r)
   *   \brief Return the r of OneStepIntegrator
   *   \param double : the value to set r
   */
  inline void setR(const int& newR)
  {
    r = newR;
  }

  /** \fn Adams* convert (OneStepIntegrator* osi)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepIntegrator* : the integrator which must be converted
   * \return a pointer on the integrator if it is of the right type, NULL otherwise
   */
  static Adams* convert(OneStepIntegrator* osi);

private:
  /** \fn Adams()
   *  \brief default constructor
   */
  Adams();

  /**  */
  int r;

};

#endif // ADAMS_H
