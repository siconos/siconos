//$Id: OneStepIntegrator.cpp,v 1.28 2005/03/08 14:23:44 jbarbier Exp $
#include "OneStepIntegrator.h"
//#include "Strategy.h"

#include "check.h"

OneStepIntegrator::OneStepIntegrator()
{
  this->integratorxml = NULL;
  this->timeDiscretisation = NULL;
  this->ds = NULL;
  this->r = 1;
}

OneStepIntegrator::OneStepIntegrator(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds)
{
  this->integratorxml = osixml;
  this->timeDiscretisation = td;
  this->ds = ds;
  cout << "  - the DynamicalSystem Linked to this OneStepIntegrator has the number " << this->ds->getNumber() << " and his id is " << this->ds->getId() << endl;
}

OneStepIntegrator::~OneStepIntegrator()
{}


void OneStepIntegrator::initialize()
{
  IN("OneStepIntegrator::initialize \n");
  this->ds->initMemory(this->r);
  OUT("OneStepIntegrator::initialize\n");
}



void OneStepIntegrator::integrate()
{
  IN("OneStepIntegrator::integrate\n");
  OUT("OneStepIntegrator::integrate\n");
}

void OneStepIntegrator::nextStep(void)
{
  IN("OneStepIntegrator::nextStep\n");
  this->ds->swapInMemory();
  OUT("OneStepIntegrator::nextStep\n");

}


void OneStepIntegrator::computeFreeState()
{
  // to do
}


void OneStepIntegrator::updateState()
{
  // to do
}


void OneStepIntegrator::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the OneStepIntegrator " << endl;
  cout << "| integratorType : " << this->integratorType << endl;
  cout << "| ds : " << this->ds->getId() << endl;
  if (this->integratorType != MOREAU_INTEGRATOR)
    cout << "| r : " << this->r << endl;
  this->timeDiscretisation->display();
  cout << "-----------------------------------------------------" << endl << endl;
}



void OneStepIntegrator::fillIntegratorWithIntegratorXML()
{
  IN("OneStepIntegrator::fillIntegratorWithIntegratorXML\n");
  if (this->integratorxml != NULL)
  {
    if (this->integratorxml->hasR()) this->r = this->integratorxml->getR();
    else cout << "Warning : the r value of the OneStepIntegrator is not defined, optional attribute." << endl;
  }
  else RuntimeException::selfThrow("OneStepIntegrator::fillIntegratorWithIntegratorXML - OneStepIntegratorXML object not exists");
  OUT("OneStepIntegrator::fillIntegratorWithIntegratorXML\n");

}

void OneStepIntegrator::saveIntegratorToXML()
{
  IN("OneStepIntegrator::saveIntegratorToXML\n");
  if (this->integratorxml != NULL)
  {
    vector<int> dsConcerned;
    dsConcerned.push_back(this->ds->getNumber());
    this->integratorxml->setDSConcerned(&dsConcerned);

    // r is saved only if the integrator is not a Moreau integrator !
    if (this->integratorType != MOREAU_INTEGRATOR) this->integratorxml->setR(this->r);
  }
  else RuntimeException::selfThrow("OneStepIntegrator::saveIntegratorToXML - OneStepIntegratorXML object not exists");
  OUT("OneStepIntegrator::saveIntegratorToXML\n");

}
//$Log: OneStepIntegrator.cpp,v $
//Revision 1.28  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.27  2005/02/28 16:22:33  jbarbier
//- rolling balls sample almost finished
//
//- in LCP, compute function now use all the interactions to make computations
//
//Revision 1.26  2005/02/14 09:52:22  charlety
//_ getters / setters put inline
//
//Revision 1.25  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.24  2005/01/18 10:35:16  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.23  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.22  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.21  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.20  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.19  2004/09/10 11:26:16  charlety
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
//Revision 1.18  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.17  2004/07/27 14:56:05  jbarbier
//- functions createStrategy, createTimeDiscretisation and createIntegrator done
//
//Revision 1.16  2004/06/29 10:38:39  acary
//Ajout des Tag CVS ID et Log
//Ajout de la gestion pas le constructeur XML de Theta
//