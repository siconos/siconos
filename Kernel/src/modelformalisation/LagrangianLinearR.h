//$Id: LagrangianLinearR.h,v 1.13 2005/02/11 17:35:55 charlety Exp $
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "Relation.h"
#include "LagrangianLinearRXML.h"
#include "SimpleVector.h"
#include "CompositeVector.h"

/** \class LagrangianLinearR
 *  \brief Lagrangian Linear Relation
 *  \author JB CHARLETY
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 * $Date: 2005/02/11 17:35:55 $
 * $Revision: 1.13 $
 * $Author: charlety $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LagrangianLinearR.h,v $
 *
 */
class LagrangianLinearR : public Relation
{
public:

  /** \fn LagrangianLinearR()
   *  \brief Default constructor
   */
  LagrangianLinearR();

  /** \fn LagrangianLinearR(RelationXML*)
   *  \brief constructor with XML object of the parent class Relation
   *  \param RelationXML* : the XML object corresponding
   */
  LagrangianLinearR(RelationXML*);

  /** \fn LagrangianLinearR(SiconosMatrix, SiconosVector)
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set h
   *  \param e SiconosVector to set b
   */
  LagrangianLinearR(SiconosMatrix, SimpleVector);
  ~LagrangianLinearR();

  /** \fn SiconosMatrix getH(void)
   *  \brief getter of the SiconosMatrix h
   *  \return a pointer on the SiconosMatrix h
   */
  inline SiconosMatrix getH(void) const
  {
    return this->h;
  } ;

  /** \fn SimpleVector getB(void)
   *  \brief getter of the SiconosVector b
   *  \return SimpleVector : value of b
   */
  inline /*SiconosVector*/SimpleVector getB(void) const
  {
    return this->b;
  };

  /** \fn SiconosMatrix* getHPtr(void)
   *  \brief getter of the SiconosMatrix* h
   *  \return a pointer on the SiconosMatrix* h
   */
  SiconosMatrix* getHPtr(void);

  /** \fn SiconosMatrix getHPtrRelatingToDS(int position)
   *  \brief getter of the SiconosMatrix H relating to one (the one at the given position in the dsVector of the interaction) of the 2 specific Dynamical System of of the Interaction
   *  \return SiconosMatrix : the H for a DS
   */
  SiconosMatrix getHRelatingToDS(int position);

  /** \fn SiconosVector* getBPtr(void)
   *  \brief getter of the SiconosVector* b
   *  \return a pointer on the SiconosVector b
   */
  SiconosVector* getBPtr(void);

  /** \fn void setH(SiconosMatrix)
   *  \brief setter on the SiconosMatrix h
   *  \param a SiconosMatrix to set h
   */
  inline void setH(const SiconosMatrix &H)
  {
    this->h = H;
  };

  /** \fn void setH(SimpleVector&)
   *  \brief set the vector b
   *  \param SimpleVector& : new value of b
   */
  inline void setB(/*SiconosVector*/SimpleVector& B)
  {
    this->b = B;
  };


  ////////////////////////////

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeFreeOutput(double time);

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeOutput(double time);

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(double time);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveRelationToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn void createRelation(LagrangianLinearRXML * relationXML,
            SiconosMatrix* H, SiconosVector* b,
            Interaction * interaction)
   *  \brief allows to create the Relation with an xml file, or the needed data
   *  \param LagrangianLinearRXML * : the XML object for this Relation
   *  \param SiconosMatrix* : the matrix H of this Relation
   *  \param SiconosVector* : the vector h of this Relation
   *  \exception RuntimeException
   */
  void createRelation(LagrangianLinearRXML * relationXML,
                      SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** \fn LagrangianLinearR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearR* convert(Relation *r);



protected:
  /** \fn void fillRelationWithRelationXML()
   *  \brief uses the RelationXML of the LagrangianLinearR to fill the fields of this Relation
   *  \exception RuntimeException
   */
  void fillRelationWithRelationXML();


private:
  /** a specific matrix to the LagrangianLinearR */
  SiconosMatrix h;

  /** a specific vector to the LagrangianLinearR */
  /*SiconosVector*/
  SimpleVector b;
};

#endif // LAGRANGIANLINEARRELATION_H
//$Log: LagrangianLinearR.h,v $
//Revision 1.13  2005/02/11 17:35:55  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.12  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.11  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.10  2005/01/31 16:26:19  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2004/09/30 08:35:02  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.8  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.7  2004/09/10 11:26:13  charlety
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
//Revision 1.6  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:37  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/07/29 14:25:36  jbarbier
//- $Log: LagrangianLinearR.h,v $
//- Revision 1.13  2005/02/11 17:35:55  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.12  2005/02/10 10:35:19  jbarbier
//- - new file regrouping all the const values of the model, modelingTools and numericalStrategy
//-
//- - new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//-
//- - new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//-
//- Revision 1.11  2005/02/04 14:52:44  jbarbier
//- - Rolling balls in progress (contact is detected)
//-
//- - time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//-
//- Revision 1.10  2005/01/31 16:26:19  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.9  2004/09/30 08:35:02  jbarbier
//- - fonction of the formalisation : fill..With...XML and link... are now
//- "protected" and no more "public"
//-
//- Revision 1.8  2004/09/22 11:16:28  charlety
//-
//- _ revision of Doxygen comments in modelformalisation
//-
//- Revision 1.7  2004/09/10 11:26:13  charlety
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
//- Revision 1.6  2004/09/03 14:41:41  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.5  2004/08/17 15:12:37  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//-
//- Revision 1.4  2004/08/12 11:55:14  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//- and $Id: LagrangianLinearR.h,v 1.13 2005/02/11 17:35:55 charlety Exp $ added
//
