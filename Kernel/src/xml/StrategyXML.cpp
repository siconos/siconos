#include "StrategyXML.h"

// includes to be deleted thanks to factories
#include "LsodarXML.h"
#include "AdamsXML.h"
#include "MoreauXML.h"
#include "LCPXML.h"
#include "QPXML.h"
#include "DFC_2DXML.h"


using namespace std;



//map<int, bool> convertVectorToMap(vector<int> v)
//{
//  cout<<"convertVectorToMap"<<endl;
//  map<int, bool> res;
//  for(int i=0; i<v.size(); i++)
//  {
//    cout<<".";
//    res[v[i]] = true;
//  }
//  return res;
//}

StrategyXML::StrategyXML()
{
  DSAvailabilityMap.clear();

  oneStepNSProblemXML = NULL;
  timeDiscretisationXML = NULL;
  /*
   * creation of the node of the TimeDiscretisation, OneStepIntegrators and the OneStepNSProblem
   */
  //  // TimeDiscretisation node
  //  xmlNewChild(this->strategyNode, NULL, (xmlChar*)STRATEGY_TIMEDISCRETISATION.c_str() , NULL);
  //
  //  // OneStepIntegrator definition node
  //  xmlNewChild(this->strategyNode, NULL, (xmlChar*)STRATEGY_ONESTEPINTEGRATOR_DEFINITION.c_str() , NULL);
  //
  //  // OneStepNSProblem node
  //  xmlNewChild(this->strategyNode, NULL, (xmlChar*)STRATEGY_ONESTEPNSPROBLEM.c_str() , NULL);
}

StrategyXML::~StrategyXML()
{
  if (timeDiscretisationXML != NULL)
    delete timeDiscretisationXML;

  if (oneStepNSProblemXML != NULL)
    delete oneStepNSProblemXML;

  if (oneStepIntegratorXMLVector.size() > 0)
  {
    for (unsigned int i = 0; i < oneStepIntegratorXMLVector.size(); i++)
    {
      delete oneStepIntegratorXMLVector[i];
    }
    oneStepIntegratorXMLVector.clear();
  }
}


StrategyXML::StrategyXML(xmlNode * rootStrategyNode, vector<int> definedNumberDS, vector<int> definedNumberInteraction)
{
  strategyNode = rootStrategyNode;
  DSAvailabilityMap.clear();
  // Fill map of available DS
  for (unsigned int i = 0; i < (definedNumberDS.size()); i++)
    DSAvailabilityMap[definedNumberDS[i]] = true; //available

  definedNumberInteractionVector = definedNumberInteraction;

  loadStrategyXML();
}


void StrategyXML::loadStrategyXML()
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, LMGC90_STRATEGY_TAG)) == NULL)
  {
    if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, ONESTEPINTEGRATOR_DEFINITION_TAG)) != NULL)
      loadOneStepIntegratorXML(node);
    else
      XMLException::selfThrow("StrategyXML - loadStrategyXML ERROR : tag " + ONESTEPINTEGRATOR_DEFINITION_TAG + " not found.");
  }
  else cout << "StrategyXML - loadStrategy : no integrators defined, use of LMGC90 tag." << endl;

  /*
   * the next node after the definition of the oneStepIntegrators is the node of the OneStepNSProblem
   */
  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, ONESTEPNSPROBLEM_TAG)) != NULL)
  {
    //if( (node = SiconosDOMTreeTools::findNodeChild(node)) != NULL )
    //    if( (SiconosDOMTreeTools::findNodeChild(node)) != NULL )
    //      {
    this->loadOneStepNSProblemXML(node);
    //      }
    //      else
    //      {
    //      XMLException::selfThrow("StrategyXML - loadStrategyXML ERROR : No OneStepNSProblem has been defined.");
    //        //cout<<"StrategyXML - loadStrategyXML ERROR : tag " << STRATEGY_ONESTEPNSPROBLEM << " not found. The OneStepNSProblem is optional"<<endl;
    //        //this->oneStepNSProblemXML = NULL;
    //      }
  }
  else
    this->oneStepNSProblemXML = NULL;

  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, TIMEDISCRETISATION_TAG)) != NULL)
    this->loadTimeDiscretisationXML(node);
  else
    XMLException::selfThrow("StrategyXML - loadStrategyXML ERROR : tag " + TIMEDISCRETISATION_TAG + " not found.");
}

void StrategyXML::loadOneStepIntegratorXML(xmlNode * rootOneStepIntegratorNode)
{
  xmlNode *node;

  string type; //Type of OneStepIntegrator
  int cpt = 0;

  if ((node = SiconosDOMTreeTools::findNodeChild((const xmlNode*)rootOneStepIntegratorNode/*, STRATEGY_ONESTEPINTEGRATOR*/)) == NULL)
    XMLException::selfThrow("StrategyXML - loadOneStepIntegratorXML ERROR : at least one " + ONESTEPINTEGRATOR_TAG + " must be declared.");

  while (node != NULL)
  {
    OneStepIntegratorXML *integratorxml;

    //type = SiconosDOMTreeTools::getStringAttributeValue(node,STRATEGY_TYPE);
    type = (char*)node->name;
    cout << "** loadOneStepIntegratorXML  " << cpt << "   type = " << type << endl;
    cpt++;
    if (type == MOREAU_TAG)
    {
      integratorxml = new MoreauXML(node, DSAvailabilityMap);
      oneStepIntegratorXMLVector.push_back(integratorxml);
    }

    else if (type == LSODAR_TAG)
    {
      integratorxml = new LsodarXML(node, DSAvailabilityMap);
      oneStepIntegratorXMLVector.push_back(integratorxml);
    }
    else if (type == ADAMS_TAG)
    {
      integratorxml = new AdamsXML(node, DSAvailabilityMap);
      oneStepIntegratorXMLVector.push_back(integratorxml);
    }
    else
      XMLException::selfThrow("StrategyXML : undefined OneStepIntegrator type : " + type);

    node = SiconosDOMTreeTools::findFollowNode(node/*, STRATEGY_ONESTEPINTEGRATOR*/);
  }
}


void StrategyXML::loadOneStepNSProblemXML(xmlNode * rootOneStepNSProblemNode)
{
  if (rootOneStepNSProblemNode == NULL)  //OneStepNSProblem is not defined
  {
    this->oneStepNSProblemXML = NULL;
  }
  else
  {
    xmlNode *node = NULL;
    /*
     * selection of the good kind of model for the OnsStepMSProblem
     * we get the first node after the OneStepNSProblem node
     * => it's the solving model node (LcpSolving, ...) ==> it's good
     * or
     * => it's the Solver node ==> the solving model node must be the following node
     */
    if (SiconosDOMTreeTools::findNodeChild(rootOneStepNSProblemNode) != NULL)
      node = SiconosDOMTreeTools::findNodeChild(rootOneStepNSProblemNode);

    //string type = SiconosDOMTreeTools::getStringAttributeValue(rootOneStepNSProblemNode, STRATEGY_TYPE);
    if (strcmp((char*)node->name, OSNSP_SOLVER.c_str()) == 0)
      node = SiconosDOMTreeTools::findFollowNode(rootOneStepNSProblemNode);

    string type((char*)node->name);
    if (type == LCP_TAG)
    {
      this->oneStepNSProblemXML = new LCPXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }
    else if (type == QP_TAG)
    {
      oneStepNSProblemXML = new QPXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }
    /*else if (type == RELAY_TAG) //--Not implemented for the moment
    {
    this->oneStepNSProblemXML= new RelayXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }*/
    else if (type == DFC_2D_TAG)
    {
      oneStepNSProblemXML = new DFC_2DXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }
    else
    {
      XMLException::selfThrow("StrategyXML : undefined OneStepNSProblem type : " + type);
    }
  }
}


void StrategyXML::loadTimeDiscretisationXML(xmlNode * rootTimeDiscretisationNode)
{
  this->timeDiscretisationXML = new TimeDiscretisationXML(rootTimeDiscretisationNode);
}

void StrategyXML::updateStrategyXML(xmlNode* node, Strategy* str)
{
  IN("StrategyXML::updateStrategyXML\n");
  this->strategyNode = node;
  this->loadStrategy(str);
  OUT("StrategyXML::updateStrategyXML\n");
}

void StrategyXML::loadStrategy(Strategy* str)
{
  IN("StrategyXML::loadStrategy\n");
  string type;
  string tmp;
  xmlNode* node;
  xmlNode* integratorDefinitionNode;
  OneStepIntegratorXML* osixml;
  OneStepNSProblemXML* osnspbxml;
  TimeDiscretisationXML* tdxml;
  int i;

  if (this->strategyNode != NULL)
  {
    // creation of the TimeDiscretisation node
    if (str->getTimeDiscretisationPtr()->getTimeDiscretisationXMLPtr() == NULL)
    {
      node = xmlNewChild(this->strategyNode, NULL, (xmlChar*)TIMEDISCRETISATION_TAG.c_str(), NULL);
      if (str->getTimeDiscretisationPtr()->isConstant())
        xmlNewProp(node, (xmlChar*)TD_ISCONSTANT.c_str(), (xmlChar*)"true");

      tdxml = new TimeDiscretisationXML();

      // linkage between the TimeDiscretisation and his TimeDiscretisationXML
      str->getTimeDiscretisationPtr()->setTimeDiscretisationXMLPtr(tdxml);

      // creation of the TimeDiscretisationXML
      tdxml->updateTimeDiscretisationXML(node, str->getTimeDiscretisationPtr());

      this->timeDiscretisationXML = tdxml;
    }

    if (SiconosDOMTreeTools::findNodeChild((const xmlNode*)this->strategyNode, LMGC90_STRATEGY_TAG) == NULL)
    {
      if (!this->hasOneStepIntegratorXML())
      {
        // creation of the integrator_Definition node
        integratorDefinitionNode = xmlNewChild(this->strategyNode, NULL, (xmlChar*)ONESTEPINTEGRATOR_DEFINITION_TAG.c_str(), NULL);
      }

      /*
       * now, creation of the OneStepIntegratorXML objects
       */
      for (i = 0; i < str->getOneStepIntegratorVectorSize(); i++)
      {
        if (str->getOneStepIntegrator(i)->getOneStepIntegratorXMLPtr() == NULL)
        {
          type = str->getOneStepIntegrator(i)->getType();
          //node = xmlNewChild( integratorDefinitionNode, NULL, (xmlChar*)STRATEGY_ONESTEPINTEGRATOR.c_str(), NULL );
          if (type == MOREAU_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)MOREAU_TAG.c_str(), NULL);
            //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_MOREAU.c_str() );
            osixml = new MoreauXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<MoreauXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            this->oneStepIntegratorXMLVector.push_back(osixml);
          }
          else if (type == LSODAR_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)LSODAR_TAG.c_str(), NULL);
            //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_LSODAR.c_str() );
            osixml = new LsodarXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<LsodarXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            this->oneStepIntegratorXMLVector.push_back(osixml);
          }
          else if (type == ADAMS_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)ADAMS_TAG.c_str(), NULL);
            //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_ADAMS.c_str() );
            osixml = new AdamsXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXMLPtr(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<AdamsXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            this->oneStepIntegratorXMLVector.push_back(osixml);
          }
          else
          {
            XMLException::selfThrow("StrategyXML - loadStrategy ERROR : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
          }
        }
      }
    }
    else
    {
      // specific treatments for the LMGC90 integrators tag
    }

    /*
     * now, creation of the OneStepNSProblemXML object if the OneStepNSProblem exists
     */
    if (str->getOneStepNSProblemPtr() != NULL)
    {
      if (!this->hasOneStepNSProblemXML())
      {
        /*
         * creation of the node "OneStepNSProblem"
         */
        node = xmlNewChild(this->strategyNode, NULL, (xmlChar*)ONESTEPNSPROBLEM_TAG.c_str(), NULL);
      }

      if (str->getOneStepNSProblemPtr()->getOneStepNSProblemXML() == NULL)
      {
        type = str->getOneStepNSProblemPtr()->getType();
        //node = xmlNewChild( this->strategyNode, NULL, (xmlChar*)STRATEGY_ONESTEPNSPROBLEM.c_str(), NULL );
        if (type == LCP_TAG)
        {
          /*node = */xmlNewChild(node, NULL, (xmlChar*)LCP_TAG.c_str(), NULL);
          //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_LCP.c_str() );
          osnspbxml = new LCPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblemPtr()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          static_cast<LCPXML*>(osnspbxml)->updateOneStepNSProblemXML(node, str->getOneStepNSProblemPtr());

          this->oneStepNSProblemXML = osnspbxml;
        }
        else if (type == QP_TAG)
        {
          /*node = */xmlNewChild(node, NULL, (xmlChar*)QP_TAG.c_str(), NULL);
          //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_LCP.c_str() );
          osnspbxml = new QPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblemPtr()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          static_cast<QPXML*>(osnspbxml)->updateOneStepNSProblemXML(node, str->getOneStepNSProblemPtr());

          this->oneStepNSProblemXML = osnspbxml;
        }
        else if (type == RELAY_TAG)
        {
          //          /*node = */xmlNewChild( node, NULL, (xmlChar*)RELAY_TAG.c_str(), NULL );
          //          osnspbxml = new RelayXML();
          //
          //          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          //          str->getOneStepNSProblem()->setOneStepNSProblemXML( osnspbxml );
          //
          //          // creation of the OneStepNSProblemXML
          //          static_cast<RelayXML*>(osnspbxml)->updateOneStepNSProblemXML( node, str->getOneStepNSProblem() );
          //
          //          this->oneStepNSProblemXML = osnspbxml;
        }
        else
        {
          XMLException::selfThrow("StrategyXML - loadStrategy ERROR : undefined OneStepNSProblem type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
        }
      }
      //      else
      //      {
      //        cout<<"StrategyXML - loadStrategy warning : redefinition of OneStepNSProblem (a OneStepNSProblem has been already defined in the XML input file but a new one has been defined in the command program)."<<endl;
      //      }
    }
  }
  else XMLException::selfThrow("StrategyXML - loadStrategy ERROR : no strategyNode defined.");
  OUT("StrategyXML::loadStrategy\n");
}

