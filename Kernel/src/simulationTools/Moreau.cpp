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
    cout << "@@@@@@@@@@@@ current time ==" << t << endl;
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

//$Log: Moreau.cpp,v $
//Revision 1.46  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.45  2005/02/28 16:22:33  jbarbier
//- rolling balls sample almost finished
//
//- in LCP, compute function now use all the interactions to make computations
//
//Revision 1.44  2005/02/24 11:03:06  charlety
//
//_ New attempt with the autotools.
//
//_ /!\ THIS VERSION IS USABLE ONLY IF YOU HAVE INSTALLED THE EXTERNAL LIBRARIES (cppunit, libxml, lapack++, lapack, nana) in /usr/ or /usr/local.
//
//_ This version was only tested on Fedora core2 for the moment.
//
//Revision 1.43  2005/02/14 09:52:21  charlety
//_ getters / setters put inline
//
//Revision 1.42  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.41  2005/01/31 16:26:26  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.40  2005/01/18 10:35:16  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.39  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.38  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.37  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.36  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.35  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.34  2004/09/10 11:26:16  charlety
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
//Revision 1.33  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.32  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.31  2004/08/12 11:55:18  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.30  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.29  2004/08/02 09:26:26  jbarbier
//- xml save for SiconosMemory corrected
//- temporary operation in Moreau::integrate because of the current version of
//SiconosVector
//
//Revision 1.28  2004/07/27 14:56:05  jbarbier
//- functions createStrategy, createTimeDiscretisation and createIntegrator done
//
//Revision 1.27  2004/07/09 14:18:55  jbarbier
//-management of the optional attributes of the DynamicalSystem
//-new node t for current time in the time of the NSDS
//
//Revision 1.26  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.25  2004/07/07 08:14:54  jbarbier
//-modifications on the test after renamming
//
//-modification of the XML schema, attributs row, col and size of Matrices and
//Vector changed from 'positiveInteger' to 'nonNegativeInteger'
//
//-function setSiconosVector/Matrix which take a SiconosVector/Matrix* in parameter to avoid
//unnecessary vector and matrix copies
//
