//$Id: LagrangianDSIO.h,v 1.4 2005/03/14 16:05:27 jbarbier Exp $
#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "DSInputOutput.h"
#include "LagrangianDSIOXML.h"

/** \class LagrangianDSIO
 *  \brief Lagrangian DSInputOutput
 *  \author Jean-Michel Barbier
 *  \version 0.1
 *  \date 17/01/2005
 *
 * $Date: 2005/03/14 16:05:27 $
 * $Revision: 1.4 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LagrangianDSIO.h,v $
 *
 */
class LagrangianDSIO : public DSInputOutput
{
public:

  /** \fn LagrangianDSIO(void);
   * \brief default constructor
   */
  LagrangianDSIO();
  ~LagrangianDSIO();

  //  /** \fn void saveDSInputOutputToXML()
  //   *  \brief copy the data of the DSInputOutput to the XML tree
  //   */
  //  void saveDSInputOutputToXML();

  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML)
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *H = NULL);

private:

  //  /** class for manage plugin (open, close librairy...) */
  //  SiconosSharedLibrary cShared;
  //
  //  /** \fn void (*computeJacobianPtr)(void);
  //   * \brief to be defined
  //   */
  //
  //  void (*computeJacobianPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* jacobPtr);
  //
  //  void (*computeHPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* yPtr);

};

#endif
//$Log: LagrangianDSIO.h,v $
//Revision 1.4  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.3  2005/03/10 12:55:20  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.2  2005/03/09 15:30:29  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.1  2005/01/17 10:56:25  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
