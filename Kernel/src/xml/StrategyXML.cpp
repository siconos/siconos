
#include "StrategyXML.h"
#include "MoreauXML.h"
#include "LsodarXML.h"
#include "AdamsXML.h"
#include "LCPXML.h"
#include "QPXML.h"
//#include "RelayXML.h"  /* not yet created */

#include "check.h"


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
  this->DSAvailabilityMap.clear();

  this->oneStepNSProblemXML = NULL;
  this->timeDiscretisationXML = NULL;
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

  if (this->oneStepIntegratorXMLVector.size() > 0)
  {
    for (int i = 0; i < this->oneStepIntegratorXMLVector.size(); i++)
    {
      delete this->oneStepIntegratorXMLVector[i];
    }
    this->oneStepIntegratorXMLVector.clear();
  }
}


StrategyXML::StrategyXML(xmlNode * rootStrategyNode, vector<int> definedNumberDS, vector<int> definedNumberInteraction)
{
  this->strategyNode = rootStrategyNode;
  this->DSAvailabilityMap.clear();
  // Fill map of available DS
  for (int i = 0; i < (definedNumberDS.size()); i++)
  {
    this->DSAvailabilityMap[definedNumberDS[i]] = true; //available
  }

  this->definedNumberInteractionVector = definedNumberInteraction;

  loadStrategyXML();
}


void StrategyXML::loadStrategyXML()
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, LMGC90_STRATEGY_TAG)) == NULL)
  {
    if ((node = SiconosDOMTreeTools::findNodeChild(strategyNode, ONESTEPINTEGRATOR_DEFINITION_TAG)) != NULL)
      this->loadOneStepIntegratorXML(node);
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
  {
    /*
       * the OneStepNSProblem is optional, it's needed only when at least an Interaction is defined
       */
    cout << "StrategyXML - loadStrategyXML Warning : tag " << ONESTEPNSPROBLEM_TAG << " not found. The OneStepNSProblem is optional" << endl;
    this->oneStepNSProblemXML = NULL;
  }

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
    cout << "** loadOneStepIntegratorXML cpt = " << cpt << "   type = " << type << endl;
    cpt++;
    if (type == MOREAU_TAG)
    {
      integratorxml = new MoreauXML(node, this->DSAvailabilityMap);
      this->oneStepIntegratorXMLVector.push_back(integratorxml);
    }

    else if (type == LSODAR_TAG)
    {
      integratorxml = new LsodarXML(node, this->DSAvailabilityMap);
      this->oneStepIntegratorXMLVector.push_back(integratorxml);
    }
    else if (type == ADAMS_TAG)
    {
      integratorxml = new AdamsXML(node, this->DSAvailabilityMap);
      this->oneStepIntegratorXMLVector.push_back(integratorxml);
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
      this->oneStepNSProblemXML = new QPXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }
    /*else if (type == RELAY_TAG) //--Not implemented for the moment
    {
      this->oneStepNSProblemXML= new RelayXML(rootOneStepNSProblemNode, definedNumberInteractionVector);
    }*/
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
    if (str->getTimeDiscretisation()->getTimeDiscretisationXML() == NULL)
    {
      node = xmlNewChild(this->strategyNode, NULL, (xmlChar*)TIMEDISCRETISATION_TAG.c_str(), NULL);
      if (str->getTimeDiscretisation()->isConstant())
        xmlNewProp(node, (xmlChar*)TD_ISCONSTANT.c_str(), (xmlChar*)"true");

      tdxml = new TimeDiscretisationXML();

      // linkage between the TimeDiscretisation and his TimeDiscretisationXML
      str->getTimeDiscretisation()->setTimeDiscretisationXML(tdxml);

      // creation of the TimeDiscretisationXML
      tdxml->updateTimeDiscretisationXML(node, str->getTimeDiscretisation());

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
        if (str->getOneStepIntegrator(i)->getOneStepIntegratorXML() == NULL)
        {
          type = str->getOneStepIntegrator(i)->getType();
          //node = xmlNewChild( integratorDefinitionNode, NULL, (xmlChar*)STRATEGY_ONESTEPINTEGRATOR.c_str(), NULL );
          if (type == MOREAU_TAG)
          {
            node = xmlNewChild(integratorDefinitionNode, NULL, (xmlChar*)MOREAU_TAG.c_str(), NULL);
            //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_MOREAU.c_str() );
            osixml = new MoreauXML();

            // linkage between the OneStepIntegrator and his OneStepIntegratorXML
            str->getOneStepIntegrator(i)->setOneStepIntegratorXML(osixml);

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
            str->getOneStepIntegrator(i)->setOneStepIntegratorXML(osixml);

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
            str->getOneStepIntegrator(i)->setOneStepIntegratorXML(osixml);

            // creation of the OneStepIntegratorXML
            static_cast<AdamsXML*>(osixml)->updateOneStepIntegratorXML(node, str->getOneStepIntegrator(i));

            this->oneStepIntegratorXMLVector.push_back(osixml);
          }
          else
          {
            XMLException::selfThrow("InteracitonXML - loadInteraction ERROR : undefined Relation type : " + type + " (have you forgotten to verify the xml files with the Siconos Schema file or update it!?).");
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
    if (str->getOneStepNSProblem() != NULL)
    {
      if (!this->hasOneStepNSProblemXML())
      {
        /*
         * creation of the node "OneStepNSProblem"
         */
        node = xmlNewChild(this->strategyNode, NULL, (xmlChar*)ONESTEPNSPROBLEM_TAG.c_str(), NULL);
      }

      if (str->getOneStepNSProblem()->getOneStepNSProblemXML() == NULL)
      {
        type = str->getOneStepNSProblem()->getType();
        //node = xmlNewChild( this->strategyNode, NULL, (xmlChar*)STRATEGY_ONESTEPNSPROBLEM.c_str(), NULL );
        if (type == LCP_TAG)
        {
          /*node = */xmlNewChild(node, NULL, (xmlChar*)LCP_TAG.c_str(), NULL);
          //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_LCP.c_str() );
          osnspbxml = new LCPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblem()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          static_cast<LCPXML*>(osnspbxml)->updateOneStepNSProblemXML(node, str->getOneStepNSProblem());

          this->oneStepNSProblemXML = osnspbxml;
        }
        else if (type == QP_TAG)
        {
          /*node = */xmlNewChild(node, NULL, (xmlChar*)QP_TAG.c_str(), NULL);
          //xmlNewProp( node, (xmlChar*)STRATEGY_TYPE.c_str(), (xmlChar*)STRATEGY_LCP.c_str() );
          osnspbxml = new QPXML();

          // linkage between the OneStepNSProblem and his OneStepNSProblemXML
          str->getOneStepNSProblem()->setOneStepNSProblemXML(osnspbxml);

          // creation of the OneStepNSProblemXML
          static_cast<QPXML*>(osnspbxml)->updateOneStepNSProblemXML(node, str->getOneStepNSProblem());

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

//$Log: StrategyXML.cpp,v $
//Revision 1.38  2005/03/23 15:03:56  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.37  2005/03/08 14:23:46  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.36  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.35  2005/01/25 09:27:18  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.34  2005/01/24 14:33:03  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.33  2005/01/13 14:14:40  jbarbier
//- correction in the XML output about size attribute in tags DS_Concerned and Interactoin _Concerned
//
//- modifications in creation of XML objects when saving data with partial XML input file
//
//Revision 1.32  2004/12/20 15:01:26  jbarbier
//- schema XML renamed V1.1
//
//- schema XML corrected about OneStepNSProblem:
//  tag OneStepNSProblem contains tags LCP, QP, ... and other tags to add
//  further
//
//Revision 1.31  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.30  2004/09/28 14:19:40  jbarbier
//
//- new test of the model : manual creation of a Model, a NSDS, some dynamical
//systems, a Strategy and some integrators.
//
//- new function : Model::saveToDOMTree(), to create/update the XML objects
//and to save the data of the platform in the DOM tree, without saving to a file.
//
//Revision 1.29  2004/09/16 11:35:26  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.28  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.27  2004/09/10 08:05:25  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.26  2004/08/20 07:34:23  jbarbier
//- creation of Model, NSDS in comand program succeed in creating SiconosModelXML,
//NSDSXML
//
//Revision 1.25  2004/08/12 14:28:38  jbarbier
//- createTimeDiscretisation in progress
//
//Revision 1.24  2004/08/05 14:35:59  charlety
//
//_ LSODAR --> Lsodar (chapter 2)
//
//Revision 1.23  2004/08/05 12:44:45  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.22  2004/07/29 14:25:46  jbarbier
