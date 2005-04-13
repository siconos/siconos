#include "Moreau.h"
#include "MoreauXML.h"
#include "LagrangianLinearTIDS.h"

#include "check.h"

Moreau::Moreau(): OneStepIntegrator()
{
  this->W = SiconosMatrix::SiconosMatrix();
  this->theta = 0.1;
  this->integratorType = MOREAU_INTEGRATOR;
}

Moreau::Moreau(OneStepIntegratorXML* osixml, TimeDiscretisation* td, DynamicalSystem* ds): OneStepIntegrator(osixml, td, ds)
{
  this->W = SiconosMatrix::SiconosMatrix();
  this->integratorType = MOREAU_INTEGRATOR;
}


Moreau::~Moreau()
{}


SiconosMatrix* Moreau::getWPtr(void)
{
  return &this->W;
}


void Moreau::initialize()
{
  IN("Moreau::initialize\n");



  OneStepIntegrator::initialize();
  this->display();

  SiconosMatrix *M;
  double h = this->timeDiscretisation->getH();

  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::initialize -- LNLDS\n"));
    RuntimeException::selfThrow("Moreau::initialize - not yet implemented for Dynamical system type :" + ds->getType());
  }
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::initialize -- LTIDS\n"));
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);
    d ->display();
    M = d->getMassPtr();

    SiconosMatrix* K = d->getKPtr();
    SiconosMatrix* C = d->getCPtr();

    // Computation of W

    this->W = *M + (h * this->theta) * (*C + (h * this->theta) * (*K));

    //\todo Taking into account of the Constraints

    // LU factorization of W with partial Pivoting
    // Warning :  The result of the factorisatoin is return in W (InPlace)

    this->W.PLUFactorizationInPlace();
  }
  else RuntimeException::selfThrow("Moreau::initialize - not yet implemented for Dynamical system type :" + ds->getType());
  // this->display();
  OUT("Moreau::initialize\n");
}


void Moreau::computeFreeState()
{
  IN("Moreau::computeFreeState\n");


  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::computeFreeState -- LNLDS\n"));
    RuntimeException::selfThrow("Moreau::computeFreeState - not yet implemented for Dynamical system type :" + ds->getType());
  }
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::computeFreeState -- LTIDS\n"));

    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);


    SiconosVector *p = d->getPPtr();
    p->zero();

    this->integrate();

    d->setQFree(d->getQ());
    d->setVelocityFree(d->getVelocity());
  }
  else RuntimeException::selfThrow("Moreau::computeFreeState - not yet implemented for Dynamical system type :" + ds->getType());

  OUT("Moreau::computeFreeState\n");
}


void Moreau::integrate()
{
  IN("Moreau::integrate()\n");

  double theta = this->theta;

  double h = this->timeDiscretisation->getH();
  double t = this->timeDiscretisation->getStrategy()->getModel()->getCurrentT();
  double told = t - h;

  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::integrate -- LNLDS\n"));
    RuntimeException::selfThrow("Moreau::integrate - not yet implemented for Dynamical system type :" + ds->getType());
  }
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::integrate -- LTIDS\n"));
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);
    SimpleVector* vold, *qold;
    SimpleVector *v, *q, *p, *x ;
    SiconosMatrix* M, *K, *C, *W;

    q = d->getQPtr();
    v = d->getVelocityPtr();

    qold = static_cast<SimpleVector*>(d->getQMemories()->getSiconosVector(0));

    vold = static_cast<SimpleVector*>(d->getVelocityMemories()->getSiconosVector(0));
    M = d->getMassPtr();
    K = d->getKPtr();
    C = d->getCPtr();
    W = this->getWPtr();
    p = d->getPPtr();

    // Inline Version
    // The method computeFExt does not allow to compute directly
    // as a function.  To do that, you have to call directly the function of the plugin
    // or call the F77 function  MoreauLTIDS

    // Computation of the Residual term
    d->computeFExt(t);
    SimpleVector FExt0 = d->getFExt();
    d->computeFExt(told);
    SimpleVector FExt1 = d->getFExt();  //

    //    *v = -h * (*C * *vold) -h * (*K * *qold) - h * h * this->theta * ( *K * *vold)
    //    + h * ((this->theta * FExt1) + (1.0 - this->theta)*FExt0) + *p;

    *v = *p;
    *v /= h;
    *v += (1.0 - this->theta) * FExt0;
    *v += (this->theta * FExt1);
    //    cout<<"Moreau::integrate ###### K then vold ################"<<endl;
    //    K->display();
    //    vold->display();
    //    cout<<"/Moreau::integrate ##################################"<<endl;
    *v -= h * this->theta * (*K * *vold);
    *v -= (*K * *qold);
    *v -= (*C * *vold);
    *v *= h;

    // Correction of the residual term to take into account Constraints if any
    // Solving of the Linear System

    this->W.PLUForwardBackwardInPlace((*v));

    *v +=  *vold;

    *q = (*qold) + h * ((this->theta * (*v)) + (1.0 - this->theta) * (*vold));
    //    *q = *v;
    //    *q *= this->theta;
    //    *q += (1.0 - this->theta) * (*vold);
    //    *q *= h;
    //    *q += *qold;


    // Right Way  : Fortran 77 version with BLAS call
    // F77NAME(MoreauLTIDS)(told,t,theta
    //                      ndof, &qold(0),&vold(0),
    //                      &W(0,0),&K(0,0),&C(0,0),fext,
    //                      &v(0),&q(0))


  }
  else RuntimeException::selfThrow("Moreau::integrate - not yet implemented for Dynamical system type :" + ds->getType());

  OUT("Moreau::integrate()\n");
}

void Moreau::updateState()
{
  IN("Moreau::updateState\n");

  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::updateState -- LNLDS\n"));
    RuntimeException::selfThrow("Moreau::updateState - not yet implemented for Dynamical system type :" + ds->getType());
  }
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::updateState -- LTIDS\n"));
    this->integrate();
  }
  else RuntimeException::selfThrow("Moreau::updateState - not yet implemented for Dynamical system type :" + ds->getType());

  OUT("Moreau::updateState\n");
}


void Moreau::display() const
{
  OneStepIntegrator::display();

  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the Moreau Integrator " << endl;
  cout << "| W " << endl;
  this->W.display();
  cout << "| theta : " << this->theta << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}


void Moreau::fillIntegratorWithIntegratorXML()
{
  IN("Moreau::fillIntegratorWithIntegratorXML\n");
  OneStepIntegrator::fillIntegratorWithIntegratorXML();
  if (this->integratorxml != NULL)
  {
    if ((static_cast<MoreauXML*>(this->integratorxml))->hasW() == true)
    {
      this->W = (static_cast<MoreauXML*>(this->integratorxml))->getW();
    }
    else
    {
      cout << "Warning : W is not defined in the XML " << endl;
    }
    if ((static_cast<MoreauXML*>(this->integratorxml))->hasTheta() == true)
    {
      this->theta = (static_cast<MoreauXML*>(this->integratorxml))->getTheta();
    }
    else
    {
      cout << "Warning :  Theta is not defined in the XML file.  Theta is set to its default value :  1.0 " << endl;
      this->theta = 1.0;
    }
  }
  else RuntimeException::selfThrow("Moreau::fillIntegratorWithIntegratorXML - IntegratorXML object not exists");
  OUT("Moreau::fillIntegratorWithIntegratorXML\n");

}

void Moreau::saveIntegratorToXML()
{
  IN("Moreau::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != NULL)
  {
    (static_cast<MoreauXML*>(this->integratorxml))->setTheta(this->theta);
    (static_cast<MoreauXML*>(this->integratorxml))->setW(&(this->W));
  }
  else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Moreau::saveIntegratorToXML\n");
}

void Moreau::saveWToXML()
{
  IN("Moreau::saveWToXML\n");
  if (this->integratorxml != NULL)
  {
    (static_cast<MoreauXML*>(this->integratorxml))->setW(&(this->W));
  }
  else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Moreau::saveWToXML\n");
}

void Moreau::createOneStepIntegrator(OneStepIntegratorXML * osiXML, TimeDiscretisation* td,
                                     DynamicalSystem* ds, double theta)//, Strategy * strategy)
{
  if (osiXML != NULL)
  {
    this->W = SiconosMatrix::SiconosMatrix();
    this->integratorxml = osiXML;
    this->timeDiscretisation = td;
    this->ds = ds;

    this->fillIntegratorWithIntegratorXML();
  }
  else
  {
    this->integratorType = MOREAU_INTEGRATOR;
    this->W = SiconosMatrix::SiconosMatrix();
    this->theta = theta;

    this->integratorxml = NULL;
    this->timeDiscretisation = td;
    this->ds = ds;
  }
}


Moreau* Moreau::convert(OneStepIntegrator* osi)
{
  cout << "Moreau::convert (OneStepIntegrator* osi)" << endl;
  Moreau* moreau = dynamic_cast<Moreau*>(osi);
  return moreau;
}
