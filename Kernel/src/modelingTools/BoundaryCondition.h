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



using namespace std;


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
//$Log: BoundaryCondition.h,v $
//Revision 1.19  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.18  2005/02/11 17:35:54  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.17  2005/01/18 17:07:37  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.16  2004/09/30 08:35:01  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.15  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.14  2004/09/10 11:26:07  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.13  2004/09/10 08:04:46  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.12  2004/08/17 15:12:37  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.11  2004/07/29 14:25:34  jbarbier
