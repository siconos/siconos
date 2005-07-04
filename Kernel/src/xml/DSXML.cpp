#include "DSXML.h"

// to be deleted thanks to factories:
#include "NLinearBCXML.h"
#include "LinearBCXML.h"
#include "PeriodicBCXML.h"
#include "LinearDSIOXML.h"
#include "LagrangianDSIOXML.h"

using namespace std;



DSXML::DSXML():
  rootDSXMLNode(NULL), parentNode(NULL), boundaryConditionXML(NULL), xMemoryXML(NULL), xDotMemoryXML(NULL), rMemoryXML(NULL),
  idNode(NULL), nNode(NULL), x0Node(NULL), xNode(NULL), xDotNode(NULL), xMemoryNode(NULL), xDotMemoryNode(NULL),
  stepsInMemoryNode(NULL), vectorFieldNode(NULL), computeJacobianXNode(NULL), boundaryConditionNode(NULL),
  dsInputOutputNode(NULL),  rNode(NULL), rMemoryNode(NULL)
{}


DSXML::DSXML(xmlNode * DSNode, const bool& isBVP):
  rootDSXMLNode(DSNode), parentNode(NULL), boundaryConditionXML(NULL), xMemoryXML(NULL), xDotMemoryXML(NULL), rMemoryXML(NULL),
  idNode(NULL), nNode(NULL), x0Node(NULL), xNode(NULL), xDotNode(NULL), xMemoryNode(NULL), xDotMemoryNode(NULL),
  stepsInMemoryNode(NULL), vectorFieldNode(NULL), computeJacobianXNode(NULL), boundaryConditionNode(NULL),
  dsInputOutputNode(NULL),  rNode(NULL), rMemoryNode(NULL)
{
  loadDSProperties(isBVP);
}

DSXML::~DSXML()
{
  if (boundaryConditionXML != NULL) delete boundaryConditionXML;
  if (rMemoryXML != NULL) delete rMemoryXML;
  if (xMemoryXML != NULL) delete xMemoryXML;
  if (xDotMemoryXML != NULL) delete xDotMemoryXML;
}

void DSXML::loadDSProperties(const bool& isBVP)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, BOUNDARYCONDITION_TAG)) != NULL)
  {
    loadBoundaryConditionXML(node);
    boundaryConditionNode = node;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, ID_ATTRIBUTE)) != NULL)
    idNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_N)) != NULL)
    nNode = node;
  else
  {
    if ((getType() != LAGRANGIAN_NON_LINEARDS_TAG) && (getType() != LAGRANGIAN_TIME_INVARIANTDS_TAG))
      XMLException::selfThrow("DSXML - loadDSProperties error : tag " + DS_N + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_X0)) != NULL)
    x0Node = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_X)) != NULL)
    xNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_XDOT)) != NULL)
    xDotNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_XMEMORY)) != NULL)
  {
    xMemoryNode = node;
    xMemoryXML = new SiconosMemoryXML(xMemoryNode, parentNode, DS_XMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_XDOTMEMORY)) != NULL)
  {
    xDotMemoryNode = node;
    xDotMemoryXML = new SiconosMemoryXML(xDotMemoryNode, parentNode, DS_XDOTMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_STEPSINMEMORY)) != NULL)
    stepsInMemoryNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_R)) != NULL)
    rNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_RMEMORY)) != NULL)
  {
    rMemoryNode = node;
    rMemoryXML = new SiconosMemoryXML(rMemoryNode, parentNode, DS_RMEMORY);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_VECTORFIELD)) != NULL)
    vectorFieldNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootDSXMLNode, DS_COMPUTEJACOBIANX)) != NULL)
    computeJacobianXNode = node;
}

void DSXML::loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode)
{
  if (rootBoundaryConditionNode == NULL)  //BoundaryCondition is not defined
  {
    boundaryConditionXML = NULL;
  }
  else
  {
    //string type = SiconosDOMTreeTools::getStringAttributeValue(rootBoundaryConditionNode, TYPE_ATTRIBUTE);
    xmlNode *node = SiconosDOMTreeTools::findNodeChild(rootBoundaryConditionNode);
    string type((char*)node->name);
    if (type == NON_LINEARBC_TAG)
    {
      boundaryConditionXML = new NLinearBCXML(node);
    }
    else if (type == LINEARBC_TAG)
    {
      boundaryConditionXML = new LinearBCXML(node);
    }
    else if (type == PERIODICBC_TAG)
    {
      boundaryConditionXML = new PeriodicBCXML(node);
    }
    else
    {
      XMLException::selfThrow("DSXML : undefined boundary condition type : " + type);
    }
  }
}

void DSXML::updateDynamicalSystemXML(xmlNode* newRootDSXMLNode, DynamicalSystem* ds, BoundaryCondition* bc)
{
  IN("DSXML::updateDynamicalSystemXML\n");
  rootDSXMLNode = newRootDSXMLNode;
  loadDS(ds);
  OUT("DSXML::updateDynamicalSystemXML\n");
}

void DSXML::loadDS(DynamicalSystem* ds)
{
  IN("DSXML::loadDS( DynamicalSystem* ds)\n");
  string type;
  xmlNode* node;

  if (ds->getBoundaryConditionPtr() != NULL)
  {
    type = ds->getBoundaryConditionPtr()->getType();
    node = xmlNewChild(rootDSXMLNode, NULL, (xmlChar*)BOUNDARYCONDITION_TAG.c_str(), NULL);
    if (type == NLINEARBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)NON_LINEARBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)NON_LINEARBC_TAG.c_str(), NULL);
      boundaryConditionXML = new NLinearBCXML();

      // linkage between the DynamicalSystem and his DSXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<NLinearBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node);  //, ds->getBoundaryCondition() );
    }
    else if (type == LINEARBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)LINEARBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)LINEARBC_TAG.c_str(), NULL);
      boundaryConditionXML = new LinearBCXML();

      // linkage between the DynamicalSystem and his DSXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<LinearBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node); //, ds->getBoundaryCondition() );
    }
    else if (type == PERIODICBC)
    {
      //xmlNewProp( node, (xmlChar*)TYPE_ATTRIBUTE.c_str(), (xmlChar*)PERIODICBC_TAG.c_str() );
      node = xmlNewChild(node, NULL, (xmlChar*)PERIODICBC_TAG.c_str(), NULL);
      boundaryConditionXML = new PeriodicBCXML();

      // linkage between the DynamicalSystem and his DSXML
      ds->getBoundaryConditionPtr()->setBoundaryConditionXML(boundaryConditionXML);

      // creation of the DynamicalSystemXML
      static_cast<PeriodicBCXML*>(boundaryConditionXML)->updateBoundaryConditionXML(node); //, ds->getBoundaryCondition() );
    }
    else
    {
      XMLException::selfThrow("DSXML - loadDS error : undefined DS type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
    }
  }

  if (ds->getDSInputOutputs().size() > 0)
  {
    int number;
    char num[32];
    map<int, DSInputOutputXML*>::iterator itdsio;
    xmlNode *dsioDefinitionNode, *nsdsNode;
    DSInputOutputXML *dsioXML;

    nsdsNode = ds->getNSDSPtr()->getNSDSXMLPtr()->getNSDSXMLNode();
    dsioDefinitionNode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)nsdsNode, DSINPUTOUTPUT_DEFINITION_TAG);
    if (dsioDefinitionNode == NULL)
      dsioDefinitionNode = xmlNewChild(nsdsNode, NULL, (xmlChar*)DSINPUTOUTPUT_DEFINITION_TAG.c_str(), NULL);

    for (unsigned int i = 0; i < ds->getDSInputOutputs().size(); i++)
    {
      if (ds->getDSInputOutput(i)->getDSInputOutputXML() == NULL)
      {
        number = ds->getDSInputOutput(i)->getNumber();
        sprintf(num, "%d", number);
        definedDSInputOutputNumbers.push_back(number);

        // verifies if this DSInputOutput has a number which not used
        itdsio = dsInputOutputXMLMap.find(number);
        if (itdsio == dsInputOutputXMLMap.end())
        {
          if (ds->getDSInputOutput(i)->getType() == LINEARDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)LINEAR_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new LinearDSIOXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else if (ds->getDSInputOutput(i)->getType() == NLINEARDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)NON_LINEAR_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new DSInputOutputXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else if (ds->getDSInputOutput(i)->getType() == LAGRANGIANDSIO)
          {
            node = xmlNewChild(dsioDefinitionNode, NULL, (xmlChar*)LAGRANGIAN_DSIO_TAG.c_str(), NULL);
            xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
            //xmlNewChild( node, NULL, (xmlChar*)DSINPUTOUTPUT_H.c_str(), NULL );
            dsioXML = new LagrangianDSIOXML();

            // linkage between the DSInputOutput and his DSInputOutputXML
            ds->getDSInputOutput(i)->setDSInputOutputXML(dsioXML);
            dsioXML->updateDSInputOutputXML(node, ds->getDSInputOutput(i));
            dsInputOutputXMLMap[number] = dsioXML;
          }
          else XMLException::selfThrow("DSXML - loadDS | Error : the DSInputOutput type : " + ds->getDSInputOutput(i)->getType() + " doesn't exist!");

          /*  end of the save : saving the DynamicalSystem linked to this DSInputOutput */
          node = xmlNewChild(node, NULL, (xmlChar*)DS_CONCERNED.c_str(), NULL);
          node = xmlNewChild(node, NULL, (xmlChar*) DYNAMICAL_SYSTEM_TAG.c_str(), NULL);
          number = ds->getNumber();
          sprintf(num, "%d", number);
          xmlNewProp(node, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
        }
        else cout << "DSXML - loadDS : the DSInputOutput type : " << ds->getDSInputOutput(i)->getType() << " already exists!" << endl;
      }
      else cout << "### strange, DSIOXML != NULL :gratgrat:" << endl;
    }
  }
  //    else
  //    cout<<"DSXML - loadDS WARNING : tag "<<DSINPUTOUTPUT_DEFINITION_TAG<<" not found,\nDefining DS InputOutput is optional."<<endl;
  OUT("DSXML::loadDS( DynamicalSystem* ds)\n");
}

DSInputOutputXML* DSXML::getDSInputOutputXML(int number)
{
  map<int, DSInputOutputXML*>::iterator it;

  it = dsInputOutputXMLMap.find(number);
  if (it == dsInputOutputXMLMap.end())
  {
    return NULL;
  }
  return dsInputOutputXMLMap[number];
}

void DSXML::setDSInputOutputXML(map<int, DSInputOutputXML*> m)
{
  definedDSInputOutputNumbers.clear();

  map<int, DSInputOutputXML*>::iterator iter;
  for (iter = m.begin(); iter != m.end(); iter++)
    definedDSInputOutputNumbers.push_back((*iter).first);

  dsInputOutputXMLMap = m;
}

//void DSXML::loadDSInputOutputXML(xmlNode * rootdsioNode)
//{
//  xmlNode *node;
//    int number; //Number of an EqualityCopnstraint
//  map<int, DSInputOutputXML*>::iterator i;
//
//    node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootdsioNode);
//
//  while(node!=NULL)
//    {
//      DSInputOutputXML *ecxml;
//
//    number = SiconosDOMTreeTools::getIntegerAttributeValue(node, NUMBER_ATTRIBUTE);
//
//    cout<<"CHARGEMENT DSInputOutput number"<<number<<endl;
//
//    i = dsInputOutputXMLMap.find(number);
//    if (i == dsInputOutputXMLMap.end())
//    {
//      definedDSInputOutputNumbers.push_back(number);
//      ecxml = new DSInputOutputXML((xmlNode *)node/*, definedDSNumbers*/);
//      dsInputOutputXMLMap[number] = ecxml;
//      }
//      else
//      {
//      XMLException::selfThrow("DSXML - loadDSInputOutputXML error : wrong DSINPUTOUTPUT number : already exists.");
//    }
//
//      node = SiconosDOMTreeTools::findFollowNode(node);
//   }
//}
