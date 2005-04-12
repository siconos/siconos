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

