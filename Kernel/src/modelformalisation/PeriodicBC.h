#ifndef PERIODICBC_H
#define PERIODICBC_H

#include "BoundaryCondition.h"

/** \class PeriodicBC
 *  \brief kind of BoundaryCondition
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) May 6, 2004
 *
 *
 *
 */

class PeriodicBC : public BoundaryCondition
{
public:

  /** \fn PeriodicBC();
   *  \brief Basic constructor
   */
  PeriodicBC();

  /** \fn PeriodicBC(BoundaryConditionXML*);
   *  \brief constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  PeriodicBC(BoundaryConditionXML*);

  ~PeriodicBC();

  /////////////////

  /** \fn void saveBCToXML()
   *  \brief copy the data of the BoundaryCondition to the XML tree
   *  \exception RuntimeException
   */
  void saveBCToXML();

  /** \fn void createBoundaryCondition(BoundaryConditionXML * bcXML)
   *  \brief allows to create the BoundaryCondition with an xml file, or the needed data
   *  \param BoundaryConditionXML* : the XML object for this BoundaryCondition
   *  \exception RuntimeException
   */
  void createBoundaryCondition(BoundaryConditionXML * bcXML);//, DynamicalSystem* ds=NULL);

  /** \fn PeriodicBC* convert (BoundaryCondition* bc)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param BoundaryCondition* : the boundary condition which must be converted
   * \return a pointer on the boundary condition if it is of the right type, NULL otherwise
   */
  static PeriodicBC* convert(BoundaryCondition* bc);

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   *  \exception RuntimeException
   */
  void fillBCWithBCXML();
};

#endif // PERIODICBC_H

//$Log: PeriodicBC.h,v $
//Revision 1.13  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.12  2005/01/31 16:26:24  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.11  2004/09/30 08:35:03  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.10  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.9  2004/09/03 14:41:47  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.8  2004/08/17 15:12:43  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.7  2004/07/29 14:25:38  jbarbier
