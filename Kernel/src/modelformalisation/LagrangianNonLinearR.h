#ifndef LAGRANGIANNONLINEARRELATION_H
#define LAGRANGIANNONLINEARRELATION_H

#include "Relation.h"
#include "LagrangianNonLinearRXML.h"

/** \class LagrangianNonLinearR
 *  \brief Lagrangian Non Linear Relation
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date Apr 27, 2004
 *
 *
 */
class LagrangianNonLinearR : public Relation
{
public:

  /** \fn LagrangianNonLinearR(void);
   * \brief default constructor
   */
  LagrangianNonLinearR();

  ~LagrangianNonLinearR();


  /** \fn void computeJacobian(void);
   * \brief default function to compute Jacobian
   */
  void computeJacobian(void);


  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  void saveRelationToXML();

  /** \fn void createRelation(LagrangianNonLinearRXML * relationXML)
   *  \brief allows to create the Relation with an xml file, or the needed data
   *  \param LagrangianNonLinearRXML * : the XML object for this Relation
   *  \param string : the name of the plugin for computeInput
   *  \param string : the name of the plugin for computeOutput
   *  \exception RuntimeException
   */
  void createRelation(LagrangianNonLinearRXML * relationXML,
                      string computeInput = "BasicPlugin:computeInput",
                      string computeOutput = "BasicPlugin:computeOutput"); //, Interaction * interaction = NULL);

  /** \fn LagrangianNonLinearR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianNonLinearR* convert(Relation *r);


private:

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** \fn void (*computeJacobianPtr)(void);
   * \brief to be defined
   */

  void (*computeJacobianPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* jacobPtr);

  void (*computeHPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* yPtr);

};

#endif // LAGRANGIANNONLINEARRELATION_H
//$Log: LagrangianNonLinearR.h,v $
//Revision 1.12  2005/02/11 17:36:00  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.11  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.10  2005/01/18 17:07:39  charlety
//
//_ added autotools makefiles for sample directory
//
//Revision 1.9  2004/12/08 12:49:37  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.8  2004/09/23 14:45:06  charlety
//
//_ Added a header file to main_siconos.cpp
//_ modified plugin functions signatures in model formalisation
//
//Revision 1.7  2004/09/23 14:09:23  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.6  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.5  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaction
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
