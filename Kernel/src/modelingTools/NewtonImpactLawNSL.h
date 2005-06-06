#ifndef NEWTONIMPACTLAWNSLAW_H
#define NEWTONIMPACTLAWNSLAW_H

#include "NonSmoothLaw.h"
#include "NewtonImpactLawNSLXML.h"

/** \class NewtonImpactLawNSL
 *  \brief Specific NonSmoothLaw for the Newton impact model
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) June 29, 2004
 *
 *
 * The class formalizes the Newton Impact law together with a complementarity condition. i.e.
 * \f[
 * \left\{\begin{array}{l}
 * y \geq 0, \lambda \geq 0, y^{T} \lambda=0\\
 *  if y \leq 0 \quad \mbox{then} \quad \dot y(t^{+}) - e \dot y(t^{-}) \geq 0,   \lambda \geq 0, (\dot y(t^{+}) - e \dot y(t^{-}))^{T} \lambda=0
 * \end{array}\right.
 * \f]
 *
 *
 */


class NewtonImpactLawNSL : public NonSmoothLaw
{

public:

  /** \fn NewtonImpactLawNSL()
   *  \brief default constructor
   */
  NewtonImpactLawNSL();

  /** \fn NewtonImpactLawNSL(NonSmoothLawXML*)
   *  \brief constructor with XML object of the NewtonImpactLawNSL
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  NewtonImpactLawNSL(NonSmoothLawXML*);

  /** \fn NewtonImpactLawNSL(const double &e)
   *  \brief constructor with the value of the NewtonImpactLawNSL attributes
   *  \param a double value e
   */
  NewtonImpactLawNSL(const double&);

  ~NewtonImpactLawNSL();

  /** \fn bool isVerified(void);
   *  \brief check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  /** \fn const double getE(void) const
   *  \brief getter of e
   *  \return the value of e
   */
  inline const double getE() const
  {
    return e;
  };

  /** \fn void setE(double)
   *  \brief setter of e
   *  \param a double to set e
   */
  inline void setE(const double& newVal)
  {
    e = newVal;
  };

  //////////////////////

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw in the XML tree
   */
  void saveNonSmoothLawToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn NewtonImpactLawNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static NewtonImpactLawNSL* convert(NonSmoothLaw* nsl);

private:
  /**  \brief The Newton coefficient of restitution
   */
  double e;
};

#endif // NewtonImpactLawNSL_H
