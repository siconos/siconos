#ifndef NSLAW_H
#define NSLAW_H

#include "NonSmoothLawXML.h"
#include "Interaction.h"

#include "SiconosConst.h"

/** \class NonSmoothLaw
 *  \brief this class contains non-smooth law
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) May 05, 2004
 *
 *
 *
 */

class Interaction;
class NonSmoothLawXML;

class NonSmoothLaw
{
public:

  /** \fn NonSmoothLaw()
   *  \brief default constructor
   */
  NonSmoothLaw();

  /** \fn NonSmoothLaw(NonSmoothLawXML*)
   *  \brief constructor with XML object of the NonSmoothLaw
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  NonSmoothLaw(NonSmoothLawXML*);

  virtual ~NonSmoothLaw();

  /** \fn bool isVerified()
   *  \brief check if the NS law is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  virtual bool isVerified(void) const = 0;

  /** \fn inline NonSmoothLawXML* getNonSmoothLawXML()
   *  \brief allows to get the NonSmoothLawXML of the NonSmoothLaw
   *  \return the pointer on the NonSmoothLawXML of the NonSmoothLaw
   */
  inline NonSmoothLawXML* getNonSmoothLawXML()
  {
    return this->nslawxml;
  }

  /** \fn inline void setNonSmoothLawXML( NonSmoothLawXML* nslawxml )
   *  \brief allows to set the NonSmoothLawXML of the NonSmoothLaw
   *  \param NonSmoothLawXML* : the pointer to set nslawxml
   */
  inline void setNonSmoothLawXML(NonSmoothLawXML* nslawxml)
  {
    this->nslawxml = nslawxml;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the NonSmoothLaw
   *  \return string : the type of the NonSmoothLaw
   */
  inline string getType() const
  {
    return this->nsLawType;
  }

  /////////////////////

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveNonSmoothLawToXML();

  /** \fn void display()
   *  \brief display the data of the NonSmoothLaw on the standard output
   *  \exception RuntimeException
   */
  virtual void display() const;


protected:
  /** \fn void fillNonSmoothLawWithNonSmoothLawXML()
   *  \brief uses the NonSmoothLawXML of the NonSmoothLaw to fill the fields of this NonSmoothLaw
   *  \exception RuntimeException
   */
  virtual void fillNonSmoothLawWithNonSmoothLawXML();


  /** type of the NonSmoothLaw */
  string nsLawType;

  /** the XML pbject linked to the NonSmoothLaw to read XML data */
  NonSmoothLawXML *nslawxml;

};

#endif // NSLAW_H
