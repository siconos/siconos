//$Id: RelayNSL.h,v 1.11 2005/02/11 17:36:02 charlety Exp $
#ifndef RELAYNSLAW_H
#define RELAYNSLAW_H

#include "NonSmoothLaw.h"
#include "RelayNSLXML.h"

/** \class RelayNSL
 *  \brief kind of Non-smooth law
 *  \author JB CHARLETY
 *  \version 0.1
 *  \date (Creation) Apr 27, 2004
 *
 * $Date: 2005/02/11 17:36:02 $
 * $Revision: 1.11 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/RelayNSL.h,v $
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
//$Log: RelayNSL.h,v $
//Revision 1.11  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.10  2005/01/31 16:26:24  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2004/09/30 08:35:03  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.8  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.7  2004/09/03 14:41:47  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.6  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.5  2004/08/12 11:55:18  jbarbier
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
//Revision 1.3  2004/07/29 14:25:38  jbarbier
//- $Log: RelayNSL.h,v $
//- Revision 1.11  2005/02/11 17:36:02  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.10  2005/01/31 16:26:24  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.9  2004/09/30 08:35:03  jbarbier
//- - fonction of the formalisation : fill..With...XML and link... are now
//- "protected" and no more "public"
//-
//- Revision 1.8  2004/09/22 11:16:28  charlety
//-
//- _ revision of Doxygen comments in modelformalisation
//-
//- Revision 1.7  2004/09/03 14:41:47  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.6  2004/08/17 15:12:43  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//-
//- Revision 1.5  2004/08/12 11:55:18  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//-
//- Revision 1.4  2004/08/11 14:43:45  jbarbier
//- - beginning of the mechanism of creation without XML input file of the objects of the platform with the
//- creatObjtect methods
//-
//- - function saveWToXML for Moreau integrator, and same specific functions to save
//- M,q and Q,p for LCP and QP
//-
//- - function to check coherency of the Model
//- and $Id: RelayNSL.h,v 1.11 2005/02/11 17:36:02 charlety Exp $ added
//
