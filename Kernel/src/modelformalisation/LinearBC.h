//$Id: LinearBC.h,v 1.18 2005/02/11 17:36:01 charlety Exp $
#ifndef LINEARBC_H
#define LINEARBC_H

#include "BoundaryCondition.h"
#include "LinearBCXML.h"

/** \class LinearBC
 *  \brief kind of BoundaryCondition
 *  \author Jean-Michel Barbier
 *  \version 1.0
 *  \date (Creation) May 6, 2004
 *
 * $Date: 2005/02/11 17:36:01 $
 * $Revision: 1.18 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LinearBC.h,v $
 *
 */
class LinearBC : public BoundaryCondition
{
public:

  /** \fn LinearBC();
   *  \brief Basic constructor
   */
  LinearBC();

  /** \fn LinearBC(BoundaryConditionXML*);
   *  \brief constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  LinearBC(BoundaryConditionXML*);

  ~LinearBC();

  /** \fn SiconosMatrix getOmega0(void)
   *  \brief allow to get the SiconosMatrix omega0
   *  \return the SiconosMatrix omega0
   */
  inline SiconosMatrix getOmega0(void) const
  {
    return this->omega0;
  };

  /** \fn SiconosMatrix getOmegaT(void)
   *  \brief allow to get the SiconosMatrix omegaT
   *  \return the SiconosMatrix omegaT
   */
  inline SiconosMatrix getOmegaT(void) const
  {
    return this->omegaT;
  };

  /** \fn SimpleVector getOmega(void)
   *  \brief get vector omega
   *  \return SimpleVector : value of omega
   */
  inline /*SiconosVector*/SimpleVector getOmega(void) const
  {
    return this->omega;
  };


  /** \fn void setOmega0(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix omega0
   *  \param the SiconosMatrix to set omega0
   */
  inline void setOmega0(SiconosMatrix &M)
  {
    this->omega0 = M;
  };

  /** \fn void setOmegaT(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix omegaT
   *  \param the SiconosMatrix to set omegaT
   */
  inline void setOmegaT(SiconosMatrix &M)
  {
    this->omegaT = M;
  };

  /** \fn void setOmega(SimpleVector&)
   *  \brief set vector omega
   *  \param SimpleVector& : new value of omega
   */
  inline void setOmega(/*SiconosVector*/SimpleVector& v)
  {
    this->omega = v;
  };


  //////////////////////

  /** \fn void saveBCToXML()
   *  \brief copy the data of the BoundaryCondition to the XML tree
   *  \exception RuntimeException
   */
  void saveBCToXML();

  /** \fn void createBoundaryCondition(BoundaryConditionXML * bcXML,
        SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
   *  \brief allows to create the BoundaryCondition with an xml file, or the needed data
   *  \param BoundaryConditionXML* : the XML object for this BoundaryCondition
   *  \param SiconosVector* : the omega vector of this BoundaryCondition
   *  \param SiconosVector* : the omega0 matrix of this BoundaryCondition
   *  \param SiconosMatrix* : the omegaT matrix of this BoundaryCondition
   *  \exception RuntimeException
   */
  void createBoundaryCondition(BoundaryConditionXML * bcXML,
                               SiconosVector* omega = NULL,
                               SiconosMatrix* omega0 = NULL, SiconosMatrix* omegaT = NULL); //,DynamicalSystem* ds=NULL);

  /** \fn LinearBC* convert (BoundaryCondition* bc)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param BoundaryCondition* : the boundary condition which must be converted
   * \return a pointer on the boundary condition if it is of the right type, NULL otherwise
   */
  static LinearBC* convert(BoundaryCondition* bc);

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   *  \exception RuntimeException
   */
  void fillBCWithBCXML();


private:
  /** Initial matrix of boundary conditions */
  SiconosMatrix omega0;
  /** Current matrix of boundary conditions */
  SiconosMatrix omegaT;
  /** Vector omega of the BoundaryCondition */
  /*SiconosVector*/
  SimpleVector omega;
};

#endif // LINEARBC_H

//$Log: LinearBC.h,v $
//Revision 1.18  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.17  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.16  2005/01/18 17:07:40  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.15  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.14  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.13  2004/09/10 11:26:14  charlety
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
//Revision 1.12  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.11  2004/08/17 15:12:39  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.10  2004/07/29 14:25:36  jbarbier
//- $Log: LinearBC.h,v $
//- Revision 1.18  2005/02/11 17:36:01  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.17  2005/01/31 16:26:20  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.16  2005/01/18 17:07:40  charlety
//-
//- _ added autotools makefiles for sample directory
//-
//- Revision 1.15  2004/09/30 08:35:02  jbarbier
//- - fonction of the formalisation : fill..With...XML and link... are now
//- "protected" and no more "public"
//-
//- Revision 1.14  2004/09/22 11:16:28  charlety
//-
//- _ revision of Doxygen comments in modelformalisation
//-
//- Revision 1.13  2004/09/10 11:26:14  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.12  2004/09/03 14:41:42  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.11  2004/08/17 15:12:39  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//- and $Id: LinearBC.h,v 1.18 2005/02/11 17:36:01 charlety Exp $ added
//
