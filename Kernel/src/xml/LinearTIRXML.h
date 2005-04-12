
/** \class LinearTIRXML
*   \brief This class manages LTIR Relation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/13/2004
*
*
*
* LinearTIRXML allows to manage data of a LTIRelation DOM tree.
*/


#ifndef __LTIRelationXML__
#define __LTIRelationXML__


#include <libxml/tree.h>

#include "RelationXML.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"


//using namespace std;


const string LTIR_C = "C";
const string LTIR_D = "D";
const string LTIR_E = "E";
const string LTIR_A = "a";


class LinearTIRXML : public RelationXML
{
public:
  LinearTIRXML();

  /** \fn LinearTIRXML(xmlNode * LTIRelationNode)
  *   \brief Build a LinearTIRXML object from a DOM tree describing a Relation with LTI type
  *   \param LinearTIRXML : the LinearTIR DOM tree
  *   \exception XMLException : if a property of the LinearTI Relation lacks in the DOM tree
  */
  LinearTIRXML(xmlNode * LTIRelationNode);

  ~LinearTIRXML();

  /** \fn SiconosMatrix getC()
  *   \brief Return the C of the LTIRelationXML
  *   \return The C SiconosMatrix of the LTIRelationXML
  */
  inline SiconosMatrix getC()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->CNode);
  }

  /** \fn SiconosMatrix getD()
  *   \brief Return the D of the LTIRelationXML
  *   \return The D SiconosMatrix of the LTIRelationXML
  */
  inline SiconosMatrix getD()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->DNode);
  }

  /** \fn SiconosMatrix getE()
  *   \brief Return the E of the LTIRelationXML
  *   \return The E SiconosMatrix of the LTIRelationXML
  */
  inline SiconosMatrix getE()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->ENode);
  }

  /** \fn SimpleVector getA()
  *   \brief Return a of the LTIRelationXML
  *   \return SimpleVector : a of LTIRelationXML
  */
  inline /*SiconosVector*/SimpleVector getA()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->aNode);
  }

  /** \fn void setC(SiconosMatrix *matrix)
  *   \brief Change the C matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for C matrix
  */
  void setC(SiconosMatrix *matrix);

  /** \fn void setD(SiconosMatrix *matrix)
  *   \brief Change the D matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for D matrix
  */
  void setD(SiconosMatrix *matrix);


  /** \fn void setE(SiconosMatrix *matrix)
  *   \brief Change the E matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for E matrix
  */
  void setE(SiconosMatrix *matrix);

  /** \fn void setA(SiconosVector *vector)
  *   \brief Change the a Vector values (in xml file or external data file switch his origin position)
  *   \param SiconosVector *vector : new value of a
  */
  void setA(SiconosVector *vector);


private:


  //Nodes
  xmlNode * CNode;
  xmlNode * DNode;
  xmlNode * ENode;
  xmlNode * aNode;

};


#endif
