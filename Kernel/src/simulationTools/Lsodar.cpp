
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


