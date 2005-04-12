
/** \class LinearDSIOXML
*   \brief This class manages LinearDSIO DSInputOutput data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 17/01/2005
*
*
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


//using namespace std;


const string LINEARDSIO_A = "A";
const string LINEARDSIO_B = "B";


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

  /** \fn SiconosMatrix getA()
  *   \brief Return the A of the LinearDSIOXML
  *   \return The A SiconosMatrix of the LinearDSIOXML
  */
  inline SiconosMatrix getA()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->ANode);
  }

  /** \fn SiconosMatrix getB()
  *   \brief Return the B of the LinearDSIOXML
  *   \return The B SiconosMatrix of the LinearDSIOXML
  */
  inline SiconosMatrix getB()
  {
    return SiconosDOMTreeTools::getSiconosMatrixValue(this->BNode);
  }

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


  /** \fn void setA(SiconosMatrix *matrix)
  *   \brief Change the A matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for A matrix
  */
  void setA(SiconosMatrix *matrix);

  /** \fn void setB(SiconosMatrix *matrix)
  *   \brief Change the B matrix values (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for B matrix
  */
  void setB(SiconosMatrix *matrix);


private:
  //Nodes
  xmlNode * ANode;
  xmlNode * BNode;
  //    xmlNode * ENode;
  //    xmlNode * aNode;
};


#endif
