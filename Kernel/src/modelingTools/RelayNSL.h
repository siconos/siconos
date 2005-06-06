#ifndef RELAYNSLAW_H
#define RELAYNSLAW_H

#include "NonSmoothLaw.h"
#include "RelayNSLXML.h"

/** \class RelayNSL
 *  \brief kind of Non-smooth law
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) Apr 27, 2004
 *
 *
 *
 * \bug
 *  \warning
 */


class RelayNSL : public NonSmoothLaw
{

public:

  /** \fn RelayNSL()
   *  \brief basic constructor
   */
  RelayNSL();

  /** \fn RelayNSL(NonSmoothLawXML*)
   *  \brief constructor with XML object of the RelayNSL
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  RelayNSL(NonSmoothLawXML*);

  /** \fn RelayNSL(const double& c, const double& d)
   *  \brief constructor with the value of the RelayNSL attributes
   *  \param a double value c
   *  \param a double value d
   */
  RelayNSL(const double&, const double&);

  ~RelayNSL();

  /** \fn bool isVerified(void);
   *  \brief check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified() const;

  /** \fn double getC()
   *  \brief to get c
   *  \return the value of c
   */
  inline const double getC() const
  {
    return c;
  };

  /** \fn void setC(const double&)
   *  \brief to set c
   *  \param a double
   */
  inline void setC(const double& newVal)
  {
    c = newVal;
  };

  /** \fn const double getD()
   *  \brief to get d
   *  \return the value of d
   */
  inline const double getD() const
  {
    return d;
  };


  /** \fn void setD(const double&)
   *  \brief to set d
   *  \param a double
   */
  inline void setD(const double& newVal)
  {
    d = newVal;
  };

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw to the XML tree
   */
  void saveNonSmoothLawToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn RelayNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static RelayNSL* convert(NonSmoothLaw* nsl);

private:
  /** represent the value after the non smooth event */
  double c;

  /** represent the value before the non smooth event */
  double d;

};

#endif // RELAYNSLAW_H
