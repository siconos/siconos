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
//$Log: ComplementarityConditionNSL.h,v $
//Revision 1.11  2005/02/11 17:35:54  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.10  2005/01/31 16:26:19  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2005/01/18 17:07:38  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.8  2004/09/30 08:35:01  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.7  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.6  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.4  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.3  2004/07/29 14:25:34  jbarbier
