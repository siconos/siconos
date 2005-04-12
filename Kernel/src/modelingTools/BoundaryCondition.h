#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <string>
#include <vector>
#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "BoundaryConditionXML.h"
#include "RuntimeException.h"
#include "SiconosConst.h"

#include "check.h"



//using namespace std;


class DynamicalSystem;

/** \class BoundaryCondition
 *  \brief represents the boundary conditions for a NSDS BVP
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) May 24, 2004
 *
 *
 */
class BoundaryCondition
{
public:

  /** \fn BoundaryCondition();
   *  \brief Basic constructor
   */
  BoundaryCondition();

  /** \fn BoundaryCondition(BoundaryConditionXML*);
   *  \brief constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  BoundaryCondition(BoundaryConditionXML*);

  virtual ~BoundaryCondition();

  /** \fn inline string getType()
   *  \brief allows to get the type of the BoundaryCondition
   *  \return string : the type of the BoundaryCondition
   */
  inline string getType() const
  {
    return this->boundaryType;
  }

  /** \fn inline void setBoundaryConditionXML( BoundaryConditionXML* bcxml )
   *  \brief allows to set the BoundaryConditionXML
   *  \param BoundaryConditionXML* : the BoundaryConditionXML of the BoundaryCondition
   */
  inline void setBoundaryConditionXML(BoundaryConditionXML* bcxml)
  {
    this->bcXML = bcxml;
  }

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   */
  virtual void fillBCWithBCXML();

  /** type of condition : linear, non linear, etc. */
  string boundaryType;
  BoundaryConditionXML* bcXML;

};

#endif // BOUNDARYCONDITION_H
