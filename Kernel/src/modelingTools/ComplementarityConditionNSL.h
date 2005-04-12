#ifndef COMPLEMENTARITYCONDITIONNSLAW_H
#define COMPLEMENTARITYCONDITIONNSLAW_H

#include "NonSmoothLaw.h"
#include "ComplementarityConditionNSLXML.h"


/** \class ComplementarityConditionNSL
 *  \brief NonSmoothLaw for complementarity models
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 *
 **/

class ComplementarityConditionNSL : public NonSmoothLaw
{
public:
  /** \fn ComplementarityConditionNSL()
   *  \brief default constructor
   */
  ComplementarityConditionNSL();

  /** \fn ComplementarityConditionNSL(NonSmoothLawXML*)
   *  \brief constructor with XML object of the parent class NonSmoothLaw
   *  \param NonSmoothLawXML* : the XML object corresponding
   */
  ComplementarityConditionNSL(NonSmoothLawXML*);

  ~ComplementarityConditionNSL();

  /** \fn bool isVerified(void);
   *  \brief checks the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  ///////////////////

  /** \fn void saveNonSmoothLawToXML()
   *  \brief copy the data of the NonSmoothLaw to the XML tree
   *  \exception RuntimeException
   */
  void saveNonSmoothLawToXML();

  /** \fn void createNonSmoothLaw(ComplementaryConditionNSLXML * nslawXML)
   *  \brief allows to create the NSLaw with an xml file, or the needed data
   *  \param ComplementarityConditionNSLXML * : the XML object for this NSLaw
   *  \exception RuntimeException
   */
  void createNonSmoothLaw(ComplementarityConditionNSLXML * nslawXML);//, Interaction * interaction = NULL);


  /** \fn ComplementarityConditionNSL* convert (NonSmoothLaw* nsl)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param NonSmoothLaw* : the law which must be converted
   * \return a pointer on the law if it is of the right type, NULL otherwise
   */
  static ComplementarityConditionNSL* convert(NonSmoothLaw* nsl);


protected:
  /** \fn void fillNonSmoothLawWithNonSmoothLawXML()
   *  \brief uses the ComplementarityConditionNSLXML of the ComplementarityConditionNSL to fill the fields of this ComplementarityConditionNSL
   *  \exception RuntimeException
   */
  void fillNonSmoothLawWithNonSmoothLawXML();

};

#endif // COMPLEMENTARITYCONDITIONNSLAW_H
