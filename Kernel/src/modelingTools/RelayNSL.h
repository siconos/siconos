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

  /** \fn RelayNSL(double c, double d)
   *  \brief constructor with the value of the RelayNSL attributes
   *  \param a double value c
   *  \param a double value d
   */
  RelayNSL(double c, double d);
  ~RelayNSL();

  /** \fn bool isVerified(void);
   *  \brief check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  /** \fn double getC(void)
   *  \brief getter of c
   *  \return the value of c
   */
  inline double getC(void) const
  {
    return this->c;
  };

  /** \fn double getD(void)
   *  \brief getter of d
   *  \return the value of d
   */
  inline double getD(void) const
  {
    return this->d;
  };

  /** \fn void setC(double)
   *  \brief setter of c
   *  \param a double to set c
   */
  inline void setC(const double C)
  {
    this->c = C;
  };

  /** \fn void setD(double)
   *  \brief setter of d
   *  \param a double to set d
   */
  inline void setD(const double D)
  {
    this->d = D;
  };


  //////////////////////

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw to the XML tree
   */
  void saveNonSmoothLawToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createNonSmoothLaw(RelayNSLXML * nslawXML, double c, double d)
   *  \brief allows to create the NSLaw with an xml file, or the needed data
   *  \param RelayNSLXML * : the XML object for this NSLaw
   *  \param double : the value for c of the RelayNSL
   *  \param double : the value for d of the RelayNSL
   *  \exception RuntimeException
   */
  void createNonSmoothLaw(RelayNSLXML * nslawXML, double c = 0, double d = 0); //, Interaction * interaction = NULL);


  /** \fn RelayNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static RelayNSL* convert(NonSmoothLaw* nsl);

protected:
  /** \fn void fillNonSmoothLawWithNonSmoothLawXML()
   *  \brief uses the RelayNSLXML of the ComplementarityConditionNSL to fill the fields of this ComplementarityConditionNSL
   *  \exception RuntimeException
   */
  void fillNonSmoothLawWithNonSmoothLawXML();


private:
  /** represent the value after the non smooth event */
  double c;

  /** represent the value before the non smooth event */
  double d;

};

#endif // RELAYNSLAW_H
