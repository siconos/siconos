//$Id: Adams.cpp,v 1.16 2005/03/08 14:23:43 jbarbier Exp $

#include "Adams.h"
#include "check.h"

Adams::Adams(): OneStepIntegrator()
{
  this->integratorType = ADAMS_INTEGRATOR;
  this->r = -1;
}

Adams::Adams(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(osixml, td, ds)
{
  this->integratorType = ADAMS_INTEGRATOR;
}

Adams::~Adams()
{}

void Adams::saveIntegratorToXML()
{
  IN("Adams::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != NULL)
  {
    //(static_cast<AdamsXML*>(this->integratorxml))->setR( this->r );
  }
  else RuntimeException::selfThrow("Adams::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Adams::saveIntegratorToXML\n");
}

void Adams::createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*td, DynamicalSystem* ds)//, Strategy * strategy)
{
  if (osiXML != NULL)
  {
    this->integratorxml = osiXML;
    this->timeDiscretisation = td;
    this->ds = ds;

    this->fillIntegratorWithIntegratorXML();
  }
  else
  {
    this->integratorxml = NULL;
    this->integratorType = ADAMS_INTEGRATOR;
    this->timeDiscretisation = td;
    this->ds = ds;
  }
}

void Adams::fillIntegratorWithIntegratorXML()
{
  IN("Adams::fillIntegratorWithIntegratorXML\n");
  OneStepIntegrator::fillIntegratorWithIntegratorXML();
  if (this->integratorxml != NULL)
  {
    if ((static_cast<AdamsXML*>(this->integratorxml))->hasR() == true)
    {
      this->r = (static_cast<AdamsXML*>(this->integratorxml))->getR();
    }
    else
    {
      cout << "Warning :  r is not defined in the XML file" << endl;
    }
  }
  else RuntimeException::selfThrow("Adams::fillIntegratorWithIntegratorXML - IntegratorXML object not exists");
  OUT("Adams::fillIntegratorWithIntegratorXML\n");
}

void Adams::initialize()
{
  IN("Adams::initialize\n");
  OneStepIntegrator::initialize();
  OUT("Adams::initialize\n");
}


Adams* Adams::convert(OneStepIntegrator* osi)
{
  cout << "Adams::convert (OneStepIntegrator* osi)" << endl;
  Adams* adams = dynamic_cast<Adams*>(osi);
  return adams;
}

//$Log: Adams.cpp,v $
//Revision 1.16  2005/03/08 14:23:43  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.15  2005/01/31 16:26:24  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.14  2004/09/23 14:09:23  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.13  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.12  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.11  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.10  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.9  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.8  2004/08/09 15:00:53  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.7  2004/07/29 14:25:39  jbarbier
//- $Log: Adams.cpp,v $
//- Revision 1.16  2005/03/08 14:23:43  jbarbier
//- - modification of constant variables :
//- in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//-
//- in simualtion tools, some constants have been moved to SiconosConst.h
//-
//- Revision 1.15  2005/01/31 16:26:24  charlety
//-
//- _ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//-
//- Revision 1.14  2004/09/23 14:09:23  jbarbier
//- - modification of the integrators, the attribute r is always optional.
//-
//- - modification of the LagrangianNonLinearR. computeInput and computeOutput are
//- required.
//-
//- Revision 1.13  2004/09/16 11:35:25  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.12  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.11  2004/09/09 08:57:44  jbarbier
//- - functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//- createTimeDiscretisation of the Strategy done.
//-
//- => all functions to create manually the objects of the platform are done
//-
//- Revision 1.10  2004/08/12 11:55:18  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//-
//- Revision 1.9  2004/08/10 14:51:49  jbarbier
//- - functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//- call the function initialize() of the base class
//-
//- Revision 1.8  2004/08/09 15:00:53  jbarbier
//- - changes in the cardinality of some attributes of the DynamicalSystem,
//- OneStepIntegrator
//-
//- - modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//-
//- - corrections in the test xml files
//- and $Id: Adams.cpp,v 1.16 2005/03/08 14:23:43 jbarbier Exp $ added
//
