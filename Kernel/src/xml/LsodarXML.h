/** \class LsodarXML
 *   \brief This class manages Lsodar data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/17/2004
 *
 *
 *
 * LsodarXML allows to manage data of a Lsodar DOM tree.
 */


#ifndef __LsodarXMLDEF__
#define __LsodarXMLDEF__

#include "OneStepIntegratorXML.h"

const std::string LSODAR_R = "r";

class LsodarXML : public OneStepIntegratorXML
{
public:

  LsodarXML();

  /** \fn LsodarXML(xmlNode * LsodarNode)
   *   \brief Build a LsodarXML object from a DOM tree describing Lsodar OneStepIntegrator
   *   \param LsodarNode : the Lsodar DOM tree
   *   \param map<int, bool> definedDSNumbers : to know if DS numbers are not used by another OneStepIntegrator
   */
  LsodarXML(xmlNode * LsodarNode,  std::map<int, bool> definedDSNumbers);


private:

  //Nodes
};


#endif
