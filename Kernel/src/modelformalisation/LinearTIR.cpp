//$Id: LinearTIR.cpp,v 1.14 2005/03/08 12:41:36 jbarbier Exp $

#include "LinearTIR.h"
#include "check.h"

LinearTIR::LinearTIR(): Relation()
{
  this->relationType = LINEARTIRELATION;//"LinearTIR";
}

LinearTIR::LinearTIR(RelationXML* relxml): Relation(relxml)
{
  this->relationType = LINEARTIRELATION;//"LinearTIR";
}

LinearTIR::~LinearTIR()
{}


SiconosMatrix* LinearTIR::getCPtr(void)
{
  return &this->C;
}

SiconosMatrix* LinearTIR::getDPtr(void)
{
  return &this->D;
}

SiconosMatrix* LinearTIR::getEPtr(void)
{
  return &this->E;
}

SiconosVector* LinearTIR::getAPtr(void)
{
  return &this->a;
}


void LinearTIR::computeOutput()
{
  IN("LinearTIR::computeOutput\n");

  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  if (vDS.size() == 2)
  {
    /*
     *\WARNING to do with CompositeVector
     */

    //    SiconosVector x(*(vDS[0]->getXPtr()), false);
    //    x.add(*( vDS[1]->getXPtr()));
  }
  else if (vDS.size() == 1)
  {
    /*
     *\WARNING to do with CompositeVector
     */
    //    SiconosVector x(*(vDS[0]->getXPtr()), false);
  }
  else RuntimeException::selfThrow("The interaction doesn't contain the right number of Dynamical Systems");

  /*
   * NOT TERMINATED... SEE LAGRANGIANLINEARRELATION
   */

  RuntimeException::selfThrow("LinearTIR::computeOutput not yet implemented");


  OUT("LinearTIR::computeOutput\n");
}



void LinearTIR::fillRelationWithRelationXML()
{
  Relation::fillRelationWithRelationXML();
  OUT("LinearTIR::fillRelationWithRelationXML\n");
  if (this->relationxml != NULL)
  {
    this->C = (static_cast<LinearTIRXML*>(this->relationxml))->getC();
    this->D = (static_cast<LinearTIRXML*>(this->relationxml))->getD();
    this->E = (static_cast<LinearTIRXML*>(this->relationxml))->getE();
    this->a = (static_cast<LinearTIRXML*>(this->relationxml))->getA();
  }
}

void LinearTIR::display() const
{
  cout << "---------------------------------------------------" << endl;
  cout << "____ data of the LinearTIR " << endl;
  cout << "| C " << endl;
  this->C.display();
  cout << "| D " << endl;
  this->D.display();
  cout << "| E " << endl;
  this->E.display();
  cout << "| a " << endl;
  this->a.display();
  cout << "____________________________" << endl;
  cout << "---------------------------------------------------" << endl;
}

void LinearTIR::saveRelationToXML()
{
  Relation::saveRelationToXML();
  OUT("LinearTIR::saveRelationToXML\n");
  if (this->relationxml != NULL)
  {
    //    this->display();

    (static_cast<LinearTIRXML*>(this->relationxml))->setC(&(this->C));
    (static_cast<LinearTIRXML*>(this->relationxml))->setD(&(this->D));
    (static_cast<LinearTIRXML*>(this->relationxml))->setE(&(this->E));
    (static_cast<LinearTIRXML*>(this->relationxml))->setA(&(this->a));
  }
}

void LinearTIR::createRelation(LinearTIRXML * relationXML,
                               SiconosMatrix* C, SiconosMatrix* D,
                               SiconosMatrix* E, SiconosVector* a)//, Interaction * interaction)
{
  if (relationXML != NULL)
  {
    this->init();
    this->relationxml = relationXML;
    this->relationType = LINEARTIRELATION;//"LinearTIR";
    this->fillRelationWithRelationXML();
  }
  else
  {
    this->C = *C;
    this->D = *D;
    this->E = *E;
    this->a = *a;
  }
}


LinearTIR* LinearTIR::convert(Relation *r)
{
  cout << "LinearTIR::convert (Relation *r)" << endl;
  LinearTIR* ltir = dynamic_cast<LinearTIR*>(r);
  return ltir;
}

//$Log: LinearTIR.cpp,v $
//Revision 1.14  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.13  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.12  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.11  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.10  2005/01/31 16:26:21  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.8  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.7  2004/09/10 11:26:15  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.6  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:39  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:16  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/07/29 14:25:37  jbarbier
//- $Log: LinearTIR.cpp,v $
//- Revision 1.14  2005/03/08 12:41:36  jbarbier
//- - constant variables files modified :
//- Some constants added in SiconosConst
//-
//- all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//-
//- Revision 1.13  2005/03/07 13:17:19  jbarbier
//- - new test : Ball2D, with a ball moving in a 2D system
//-
//- - another constant variables moved/refactored in XMLTagsName
//- - making uniform the name of the constant variables
//-
//- Revision 1.12  2005/02/11 17:36:01  charlety
//-
//- _ little "inspection of code"
//- _ basic getters and setters passed inline
//- _ getters functions passed const
//-
//- Revision 1.11  2005/02/01 11:08:42  charlety
//-
//- _ some displays of values during computations suppressed.
//-
//- Revision 1.10  2005/01/31 16:26:21  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.9  2004/09/22 11:16:28  charlety
//-
//- _ revision of Doxygen comments in modelformalisation
//-
//- Revision 1.8  2004/09/16 11:35:24  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.7  2004/09/10 11:26:15  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.6  2004/09/03 14:41:42  jbarbier
//- - new functions to create the boundary condition of the dynamical systems
//- - new functions to add an interaction to a NSDS
//- - new functions to create the relation and the non-smooth law of an interaciton
//-
//- Revision 1.5  2004/08/17 15:12:39  jbarbier
//- - methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//- createRelation and createNSLaw completed with the required attributes
//-
//- Revision 1.4  2004/08/12 11:55:16  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//- and $Id: LinearTIR.cpp,v 1.14 2005/03/08 12:41:36 jbarbier Exp $ added
//
