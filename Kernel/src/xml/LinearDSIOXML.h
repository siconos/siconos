//$Id: LinearDSIOXML.h,v 1.4 2005/03/10 12:55:22 jbarbier Exp $

/** \class LinearDSIOXML
*   \brief This class manages LinearDSIO DSInputOutput data
*   \author Jean-Michel Barbier
*   \version 1.0
*   \date 17/01/2005
*
*
* $Date: 2005/03/10 12:55:22 $
* $Revision: 1.4 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/LinearDSIOXML.h,v $
*
* LinearTIRXML allows to manage data of a LinearDSIO in the DOM tree.
*/


#ifndef __LinearDSIOXML__
#define __LinearDSIOXML__


#include <libxml/tree.h>

#include "DSInputOutputXML.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"


using namespace std;


const string LINEARDSIO_C = "C";
const string LINEARDSIO_D = "D";
const string LINEARDSIO_E = "E";
const string LINEARDSIO_A = "a";


class LinearDSIOXML : public DSInputOutputXML
{
public:
  LinearDSIOXML();

  /** \fn LinearDSIOXML(xmlNode * , vector<int> )
  *   \brief Build a DSInputOutputXML object from a DOM tree describing a DSInputOutput
  *   \param xmlNode* : the DSInputOutput DOM tree
  //    *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the DSInputOutput (identified by number) exists
  */
  LinearDSIOXML(xmlNode * dsioNode/*, vector<int> definedDSNumbers */);
  ~LinearDSIOXML();

  //    /** \fn SiconosMatrix getC()
  //    *   \brief Return the C of the LinearDSIOXML
  //    *   \return The C SiconosMatrix of the LinearDSIOXML
  //    */
  //    inline SiconosMatrix getC()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->CNode);
  //    }
  //
  //    /** \fn SiconosMatrix getD()
  //    *   \brief Return the D of the LinearDSIOXML
  //    *   \return The D SiconosMatrix of the LinearDSIOXML
  //    */
  //    inline SiconosMatrix getD()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->DNode);
  //    }
  //
  //    /** \fn SiconosMatrix getE()
  //    *   \brief Return the E of the LinearDSIOXML
  //    *   \return The E SiconosMatrix of the LinearDSIOXML
  //    */
  //    inline SiconosMatrix getE()
  //    {
  //      return SiconosDOMTreeTools::getSiconosMatrixValue(this->ENode);
  //    }
  //
  //    /** \fn SimpleVector getA()
  //    *   \brief Return a of the LinearDSIOXML
  //    *   \return SimpleVector : a of LinearDSIOXML
  //    */
  //    inline /*SiconosVector*/SimpleVector getA()
  //    {
  //      return SiconosDOMTreeTools::getSiconosVectorValue(this->aNode);
  //    }
  //
  //    /** \fn void setC(SiconosMatrix *matrix)
  //    *   \brief Change the C matrix values (in xml file or external data file switch his origin position)
  //    *   \param SiconosMatrix matrix : the new value for C matrix
  //    */
  //    void setC(SiconosMatrix *matrix);
  //
  //    /** \fn void setD(SiconosMatrix *matrix)
  //    *   \brief Change the D matrix values (in xml file or external data file switch his origin position)
  //    *   \param SiconosMatrix matrix : the new value for D matrix
  //    */
  //    void setD(SiconosMatrix *matrix);
  //
  //
  //    /** \fn void setE(SiconosMatrix *matrix)
  //    *   \brief Change the E matrix values (in xml file or external data file switch his origin position)
  //    *   \param SiconosMatrix matrix : the new value for E matrix
  //    */
  //    void setE(SiconosMatrix *matrix);
  //
  //    /** \fn void setA(SiconosVector *vector)
  //    *   \brief Change the a Vector values (in xml file or external data file switch his origin position)
  //    *   \param SiconosVector *vector : new value of a
  //    */
  //    void setA(SiconosVector *vector);


private:
  //Nodes
  //    xmlNode * CNode;
  //    xmlNode * DNode;
  //    xmlNode * ENode;
  //    xmlNode * aNode;
};


#endif
//$Log: LinearDSIOXML.h,v $
//Revision 1.4  2005/03/10 12:55:22  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.3  2005/03/09 15:30:37  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.2  2005/01/26 13:50:40  jbarbier
//
//- loading of an XML input file now loads EqualityConstraints and DSInputOutputs
//
//Revision 1.1  2005/01/17 10:56:27  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
