//$Id: TimeDiscretisationXML.cpp,v 1.13 2004/09/16 11:35:26 jbarbier Exp $

#include "TimeDiscretisationXML.h"

#include "check.h"


TimeDiscretisationXML::TimeDiscretisationXML()
{
  this->hMin = false;
  this->hMax = false;
  this->hMinNode = false;
  this->hMaxNode = false;
  this->hNode = NULL;
  this->NNode = NULL;
  this->tkNode = NULL;
}

TimeDiscretisationXML::TimeDiscretisationXML(xmlNode * timeDiscretisationNode)
{
  this->loadTimeDiscretisationProperties(timeDiscretisationNode);
}


void TimeDiscretisationXML::loadTimeDiscretisationProperties(xmlNode * timeDiscretisationNode)
{
  xmlNode *node;

  rootNode = timeDiscretisationNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_H)) != NULL)
  {
    this->hNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_H + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_H << " not found.Optional attribute" << endl;
    this->hNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_N)) != NULL)
  {
    this->NNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_N + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_N << " not found.Optional attribute" << endl;
    this->NNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_TK)) != NULL)
  {
    this->tkNode = node;
  }
  else
  {
    //XMLException::selfThrow("TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " + TD_TK + " not found.");
    cout << "TimeDiscretisationXML - loadTimeDiscretisationProperties error : tag " << TD_TK << " not found.Optional attribute" << endl;
    this->tkNode = NULL;
  }

  //Optional attribute
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMIN)) != NULL)
  {
    this->hMinNode = node;
    this->hMin = true;
  }
  else this->hMin = false;
  //Optional attribute
  if ((node = SiconosDOMTreeTools::findNodeChild(timeDiscretisationNode, TD_HMAX)) != NULL)
  {
    this->hMaxNode = node;
    this->hMax = true;
  }
  else this->hMax = false;
}

void TimeDiscretisationXML::updateTimeDiscretisationXML(xmlNode* node, TimeDiscretisation* td)
{
  IN("TimeDiscretisation::updateTimeDiscretisationXML\n");
  this->rootNode = node;
  //  this->loadNSDS( nsds );
  OUT("TimeDiscretisation::updateTimeDiscretisationXML\n");
}

//$Log: TimeDiscretisationXML.cpp,v $
//Revision 1.13  2004/09/16 11:35:26  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.12  2004/09/15 13:23:14  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.11  2004/08/05 09:31:38  jbarbier
//- test successfull on the TimeDiscretisation for the triplet (T, t0, N),
//...
//
//- cjecking of this triplet in the createTimeDiscretisation method
//
//Revision 1.10  2004/07/29 14:25:49  jbarbier
//- $Log: TimeDiscretisationXML.cpp,v $
//- Revision 1.13  2004/09/16 11:35:26  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.12  2004/09/15 13:23:14  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.11  2004/08/05 09:31:38  jbarbier
//- - test successfull on the TimeDiscretisation for the triplet (T, t0, N),
//- ...
//-
//- - cjecking of this triplet in the createTimeDiscretisation method
//- and $Id: TimeDiscretisationXML.cpp,v 1.13 2004/09/16 11:35:26 jbarbier Exp $ added
//
