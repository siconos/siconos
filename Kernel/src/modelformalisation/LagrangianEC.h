//$Id: LagrangianEC.h,v 1.5 2005/03/15 09:57:47 jbarbier Exp $
#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "EqualityConstraint.h"
#include "LagrangianECXML.h"

/** \class LagrangianEC
 *  \brief Lagrangian EqualityConstraint
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/03/15 09:57:47 $
 * $Revision: 1.5 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LagrangianEC.h,v $
 *
 */
class LagrangianEC : public EqualityConstraint
{
public:

  /** \fn LagrangianEC(void);
   * \brief default constructor
   */
  LagrangianEC();
  ~LagrangianEC();

  LagrangianEC(EqualityConstraintXML*);


  /** \fn void createDSInputOutput(EqualityConstraintXML * ecXML)
   *  \brief allows to create the EqualityConstraint with an xml file, or the needed data
   *  \param LagrangianECXML * : the XML object for this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, vector<DSInputOutput*> *dsioVector = NULL);

private:

};

#endif
//$Log: LagrangianEC.h,v $
//Revision 1.5  2005/03/15 09:57:47  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.4  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.3  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.2  2005/03/10 12:55:20  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.1  2005/03/09 15:30:30  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
