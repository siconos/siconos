//$Id: NLinearBC.h,v 1.13 2005/02/11 17:36:01 charlety Exp $
#ifndef NLINEARBC_H
#define NLINEARBC_H

#include "BoundaryCondition.h"

/** \class NLinearBC
 *  \brief kind of BoundaryCondition
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date (Creation) May 6, 2004
 *
 * $Date: 2005/02/11 17:36:01 $
 * $Revision: 1.13 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/NLinearBC.h,v $
 *
 *
 */
class NLinearBC : public BoundaryCondition
{
public:

  /** \fn NLinearBC()
   *  \brief default constructor
   */
  NLinearBC();

  /** \fn NLinearBC(BoundaryConditionXML*)
   *  \brief constructor with XML object of the boundary condition
   *  \param BoundaryConditionXML* : the XML object corresponding
   */
  NLinearBC(BoundaryConditionXML*);

  ~NLinearBC();


  /////////////////////

  /** \fn void saveBCToXML()
   *  \brief copy the data of the BoundaryCondition to the XML tree
   */
  void saveBCToXML();

  /** \fn void createBoundaryCondition(BoundaryConditionXML * bcXML)
   *  \brief allows to create the BoundaryCondition with an xml file, or the needed data
   *  \param BoundaryConditionXML* : the XML object for this BoundaryCondition
   *  \exception RuntimeException
   */
  void createBoundaryCondition(BoundaryConditionXML * bcXML);

  /** \fn NLinearBC* convert (BoundaryCondition* bc)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param BoundaryCondition* : the boundary condition which must be converted
   * \return a pointer on the boundary condition if it is of the right type, NULL otherwise
   */
  static NLinearBC* convert(BoundaryCondition* bc);

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   */
  void fillBCWithBCXML();
};

#endif // NLINEARBC_H


//$Log: NLinearBC.h,v $
//Revision 1.13  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.12  2005/01/31 16:26:21  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.11  2005/01/18 17:07:41  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.10  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.9  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.8  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.7  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.6  2004/07/29 14:25:37  jbarbier
//- $Log: NLinearBC.h,v $
//- Revision 1.13  2005/02/11 17:36:01  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.12  2005/01/31 16:26:21  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.11  2005/01/18 17:07:41  charlety
//-
//- _ added autotools makefiles for sample directory
//-
//- Revision 1.10  2004/09/30 08:35:02  jbarbier
//- - fonction of the formalisation : fill..With...XML and link... are now
//- "protected" and no more "public"
//-
//- Revision 1.9  2004/09/22 11:16:28  charlety
//-
//- _ revision of Doxygen comments in modelformalisation
//-
//- Revision 1.8  2004/09/03 14:41:42  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.7  2004/08/17 15:12:43  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//- and $Id: NLinearBC.h,v 1.13 2005/02/11 17:36:01 charlety Exp $ added
//
