//$Id: LinearDSIO.h,v 1.4 2005/03/14 16:05:27 jbarbier Exp $
#ifndef LINEARDSIO_H
#define LINEARDSIO_H

#include "DSInputOutput.h"
#include "LinearDSIOXML.h"

/** \class LinearDSIO
 *  \brief Linear DSInputOutput
 *  \author JB CHARLETY
 *  \version 0.1
 *  \date Apr 27, 2004
 *
 * $Date: 2005/03/14 16:05:27 $
 * $Revision: 1.4 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelformalisation/LinearDSIO.h,v $
 *
 *
 */
class LinearDSIO : public DSInputOutput
{
public:

  /** \fn LinearDSIO();
   *  \brief Basic constructor
   */
  LinearDSIO();

  /** \fn LinearDSIO(LinearDSIOXML*)
   *  \brief constructor with XML object of the LinearDSIO
   *  \param LinearDSIOXML* : the XML object corresponding
   */
  LinearDSIO(DSInputOutputXML*);

  ~LinearDSIO();

  //  /** \fn void saveDSInputOutputToXML()
  //   *  \brief copy the data of the DSInputOutput to the XML tree
  //   */
  //  void saveDSInputOutputToXML();
  //
  /** \fn void createDSInputOutput(DSInputOutputXML * dsioXML)
   *  \brief allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *H = NULL);


protected:
  //  /** \fn void fillDSInputOutputWithDSInputOutputXML()
  //   *  \brief uses the DSInputOutputXML of the LinearDSIO to fill the fields of this DSInputOutput
  //   */
  //  void fillDSInputOutputWithDSInputOutputXML();


private:
};

#endif // LINEARDSIO_H
//$Log: LinearDSIO.h,v $
//Revision 1.4  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.3  2005/03/10 12:55:20  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.2  2005/03/09 15:30:31  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.1  2005/01/17 10:56:25  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
