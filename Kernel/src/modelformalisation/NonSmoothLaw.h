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
//$Log: NonSmoothLaw.h,v $
//Revision 1.8  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.7  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.6  2004/09/30 08:35:03  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.5  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.4  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.3  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.2  2004/07/29 14:25:37  jbarbier
