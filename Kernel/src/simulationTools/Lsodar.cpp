
#include "Lsodar.h"
#include "AdamsXML.h"

Lsodar::Lsodar(): OneStepIntegrator()
{
  this->integratorType = LSODAR_INTEGRATOR;
}

Lsodar::Lsodar(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(osixml, td, ds)
{
  this->integratorType = LSODAR_INTEGRATOR;
}

Lsodar::~Lsodar()
{}

void Lsodar::computeFreeState()
{
  IN("Lsodar::computeFreeState\n");

  this->integrate();

  OUT("Lsodar::computeFreeState\n");
}

void Lsodar::integrate()
{
  IN("Lsodar::integrate\n");

  //testFunctionInParameter(ds->getVectorFieldPtr());
  tryfunction(ds->getVectorFieldPtr());

  OUT("Lsodar::integrate\n");
}


void Lsodar::saveIntegratorToXML()
{
  IN("Lsodar::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != NULL)
  {
    //(static_cast<LsodarXML*>(this->integratorxml))->setR( this->r );
  }
  else RuntimeException::selfThrow("Adams::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Lsodar::saveIntegratorToXML\n");
}

void Lsodar::createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation*td, DynamicalSystem* ds)//, Strategy * strategy)
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
    this->integratorType = LSODAR_INTEGRATOR;
    this->timeDiscretisation = td;
    this->ds = ds;
  }
}



void tryfunction(fctPtr f)
{
  printf("in testFunctionInParameter \n") ;

  //typedef void (* fctPtr)(int *sizeOfX, double *time, double *x, double *xdot);

  //  double x = [2];
  //  double time = 0;
  //  double xDot[2] = 0;
  //  int sizeOfX = 2;

  double *x;
  double *time;
  double *xDot;
  int *sizeOfX;

  x = (double*)malloc(sizeof(double) * 2);
  xDot = (double*)malloc(sizeof(double) * 2);
  time = (double*)malloc(sizeof(double));
  sizeOfX = (int*)malloc(sizeof(int));

  *sizeOfX = 2;
  x[0] = 1;
  x[1] = -1;
  printf("sizeOfX = %i -- %d\n", *sizeOfX, *sizeOfX);

  f(sizeOfX, time, x, xDot);
  printf("out testFunctionInParameter \n")  ;

  free(x);
  free(xDot);
  free(time);
  free(sizeOfX);
}


void Lsodar::fillIntegratorWithIntegratorXML()
{
  IN("Lsodar::fillIntegratorWithIntegratorXML\n");
  OneStepIntegrator::fillIntegratorWithIntegratorXML();
  if (this->integratorxml != NULL)
  {
    if ((static_cast<LsodarXML*>(this->integratorxml))->hasR() == true)
    {
      this->r = (static_cast<AdamsXML*>(this->integratorxml))->getR();
    }
    else
    {
      cout << "Warning :  r is not defined in the XML file" << endl;
    }
  }
  else RuntimeException::selfThrow("Adams::fillIntegratorWithIntegratorXML - IntegratorXML object not exists");
  OUT("Lsodar::fillIntegratorWithIntegratorXML\n");
}

void Lsodar::initialize()
{
  IN("Lsodar::initialize\n");
  OneStepIntegrator::initialize();
  OUT("Lsodar::initialize\n");
}

Lsodar* Lsodar::convert(OneStepIntegrator* osi)
{
  cout << "Lsodar::convert (OneStepIntegrator* osi)" << endl;
  Lsodar* lsodar = dynamic_cast<Lsodar*>(osi);
  return lsodar;
}


//$Log: Lsodar.cpp,v $
//Revision 1.13  2005/03/08 14:23:43  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.12  2005/01/31 16:26:25  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.11  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.10  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.9  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.8  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.7  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.6  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.5  2004/08/10 14:51:49  jbarbier
//- functions initialize() of the Lsodar and Adams OneStepIntegrator completed to
//call the function initialize() of the base class
//
//Revision 1.4  2004/08/10 13:04:19  charlety
//
//_ try to pass a pointer of a C function of th eplugin to a numerical routine
//
//Revision 1.3  2004/08/09 15:00:54  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.2  2004/08/06 10:46:30  charlety
//
//_ example Oscillator in progress
//_ corrected a bug : theXML  save of a system without OneStepIntegrator was not OK.
//
//Revision 1.1  2004/08/05 14:33:21  charlety
//
//_ class LSODAR is now named Lsodar
//
//Revision 1.7  2004/07/29 14:25:39  jbarbier
