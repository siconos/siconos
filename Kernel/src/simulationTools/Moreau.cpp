#include "Moreau.h"
#include "MoreauXML.h"
#include "LagrangianLinearTIDS.h"

#include "check.h"

// --- xml constructor ---
Moreau::Moreau(OneStepIntegratorXML *osiXML, TimeDiscretisation* td, DynamicalSystem* ds) : OneStepIntegrator(td, ds), W(0), theta(0.1)
{
  integratorxml = osiXML;
  integratorType = MOREAU_INTEGRATOR;
  // Memory allocation for W
  int sizeW = ds->getXPtr()->size();
  cout << " TAILLE W" << sizeW << endl;
  W = new SiconosMatrix(sizeW, sizeW);

  // xml loading
  if (osiXML != 0)
  {
    if ((static_cast<MoreauXML*>(this->integratorxml))->hasW() == true)
    {
      *this->W = (static_cast<MoreauXML*>(this->integratorxml))->getW();
    }
    if ((static_cast<MoreauXML*>(this->integratorxml))->hasTheta() == true)
    {
      this->theta = (static_cast<MoreauXML*>(this->integratorxml))->getTheta();
    }
  }
  else RuntimeException::selfThrow("Moreau::Moreau() - xml constructor - IntegratorXML object not exists");
}

// --- constructor from a minimum set of data ---
Moreau::Moreau(TimeDiscretisation* td, DynamicalSystem* ds, const double& newTheta): OneStepIntegrator(td, ds), W(0), theta(newTheta)
{
  integratorType = MOREAU_INTEGRATOR;
  // Memory allocation for W
  int sizeW = ds->getXPtr()->size();
  cout << " TAILLE W" << sizeW << endl;
  W = new SiconosMatrix(sizeW, sizeW);
}

Moreau::~Moreau()
{
  delete W;
  W = 0;
}

void Moreau::initialize()
{
  IN("Moreau::initialize\n");
  OneStepIntegrator::initialize();

  // General Lagrangian system
  if (this->ds->getType() == LNLDS)
  {
    // nothing to do: W is computed at each step
  }
  // Lagrangian linear time invariant system: computation of W Moreau's iteration matrix
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::initialize -- LTIDS\n"));
    // get the ds
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);
    // get pointer on mass matrix
    SiconosMatrix *M;
    M = d->getMassPtr();
    // get K and C pointers
    SiconosMatrix* K = d->getKPtr();
    SiconosMatrix* C = d->getCPtr();
    // get time step
    double h = this->timeDiscretisation->getH();

    // --- Computation of W ---
    *W = *M + (h * this->theta) * (*C + (h * this->theta) * (*K));

    //\todo Taking into account of the Constraints

    // LU factorization of W with partial Pivoting
    // Warning :  The result of the factorisatoin is return in W (InPlace)

    W->PLUFactorizationInPlace();
  }
  else RuntimeException::selfThrow("Moreau::initialize - not yet implemented for Dynamical system type :" + ds->getType());
  OUT("Moreau::initialize\n");
}


void Moreau::computeFreeState()
{
  IN("Moreau::computeFreeState\n");
  // general Lagrangian system
  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::computeFreeState -- LNLDS\n"));
    // State i+1: present time step
    // index k corresponds to Newton known step => actual saved state in DS

    // --- Wk calculous ---
    // -- Get the DS --
    LagrangianDS* d = static_cast<LagrangianDS*>(this->ds);
    // get current time, theta and time step
    double t = this->timeDiscretisation->getStrategy()->getModel()->getCurrentT();
    double theta = this->theta;
    double h = this->timeDiscretisation->getH();

    // --- Computing of Moreau's matrix of iteration ---
    //
    // compute Mass matrix for state i+1,k
    d->computeMass(t);
    SiconosMatrix *M = d->getMassPtr();
    // compute Kt matrix for state i+1,k
    d->computeJacobianQFInt(t);
    SiconosMatrix *JacoQFInt = d->getJacobianQFIntPtr();
    d->computeJacobianQQNLInertia(t);
    SiconosMatrix *JacoQQNL = d->getJacobianQQNLInertiaPtr();
    // compute Ct matrix for state i+1,k
    d->computeJacobianVelocityFInt(t);
    SiconosMatrix *JacoVFInt = d->getJacobianVelocityFIntPtr();
    d->computeJacobianVelocityQNLInertia(t);
    SiconosMatrix *JacoVQNL = d->getJacobianVelocityQNLInertiaPtr();
    // calculate Wk
    *W = *M + h * theta * ((*JacoVFInt) + (*JacoVQNL) + h * theta * (*JacoQFInt) + (*JacoQQNL));
    // LU factorization of W
    W->PLUFactorizationInPlace();

    // --- RESfree calculus ---
    //
    // Get state i (previous time step)
    SimpleVector* vold, *qold;
    qold = static_cast<SimpleVector*>(d->getQMemoryPtr()->getSiconosVector(0));
    vold = static_cast<SimpleVector*>(d->getVelocityMemoryPtr()->getSiconosVector(0));
    // Previous time step (i)
    double told = t - h;
    // Compute QNL, Fint and Fext for state i
    // warning: get values and not pointers
    d->computeQNLInertia(qold, vold);
    SimpleVector QNL0 = d->getQNLInertia();
    d->computeFInt(told, qold, vold);
    SimpleVector FInt0 = d->getFInt();
    d->computeFExt(told);
    SimpleVector FExt0 = d->getFExt();
    // Compute QNL, Fint and Fext for present state
    // warning: get values and not pointers
    d->computeQNLInertia();
    SimpleVector QNL1 = d->getQNLInertia();
    d->computeFInt(t);
    SimpleVector FInt1 = d->getFInt();
    d->computeFExt(t);
    SimpleVector FExt1 = d->getFExt();
    // RESfree ...
    SimpleVector *v = d->getVelocityPtr();
    SimpleVector *RESfree = new SimpleVector(FExt1.size());
    *RESfree = *M * (*v - *vold) + h * ((1.0 - theta) * (QNL0 + FInt0 - FExt0) + theta * (QNL1 + FInt1 - FExt1));

    // ------- velocity free -------------
    // Solution of Wk(v-vfree)=RESfree
    SimpleVector *vfree = d->getVelocityFreePtr();
    *vfree = *v - W->PLUForwardBackward((*RESfree));

    // calculate qfree (whereas it is useless for future computation)
    SimpleVector *qfree = d->getQFreePtr();
    *qfree = (*qold) + h * (theta * (*vfree) + (1.0 - theta) * (*vold));

    delete RESfree;
    // update q and velocity of present state with free state values (and not pointers !!): useless, just to be coherent with LTIDS case
    d->setQ(d->getQFree());
    d->setVelocity(d->getVelocityFree());
  }

  // Lagrangian linear time invariant case
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::computeFreeState -- LTIDS\n"));
    // get dynamical system
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);

    // get p pointer and initialize it to zero
    SimpleVector *p = d->getPPtr();
    p->zero();

    // integration -> free state computation and save in present state
    this->integrate();

    // save results of integrate (q and velocity) in free state
    d->setQFree(d->getQ());
    d->setVelocityFree(d->getVelocity());
  }
  else RuntimeException::selfThrow("Moreau::computeFreeState - not yet implemented for Dynamical system type :" + ds->getType());
  // Remark: at the end of this step, present state and free state are equal.
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
    //VL(("Moreau::integrate -- LNLDS\n"));
    // We do not use integrate() for LNDS
  }
  else if (this->ds->getType() == LTIDS)
  {
    VL(("Moreau::integrate -- LTIDS\n"));
    // get the ds
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(this->ds);
    // get q and velocity pointers for current time step
    SimpleVector *v, *q, *vold, *qold;
    q = d->getQPtr();
    v = d->getVelocityPtr();
    // get q and velocity pointers for previous time step
    qold = static_cast<SimpleVector*>(d->getQMemoryPtr()->getSiconosVector(0));
    vold = static_cast<SimpleVector*>(d->getVelocityMemoryPtr()->getSiconosVector(0));
    // get mass, K and C pointers
    SiconosMatrix *M, *K, *C;
    M = d->getMassPtr();
    K = d->getKPtr();
    C = d->getCPtr();
    // get p pointer
    SimpleVector  *p;
    p = d->getPPtr();

    // Inline Version
    // The method computeFExt does not allow to compute directly
    // as a function.  To do that, you have to call directly the function of the plugin
    // or call the F77 function  MoreauLTIDS

    // Computation of the external forces
    d->computeFExt(told);
    SimpleVector FExt0 = d->getFExt();
    d->computeFExt(t);
    SimpleVector FExt1 = d->getFExt();

    // velocity computation
    *v = h * (this->theta * FExt1 + (1.0 - this->theta) * FExt0 - (*C * *vold) - (*K * *qold) - h * this->theta * (*K * *vold)) + *p;
    W->PLUForwardBackwardInPlace((*v));
    *v +=  *vold;

    // q computation
    *q = (*qold) + h * ((this->theta * (*v)) + (1.0 - this->theta) * (*vold));

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
  // general Lagrangian system
  if (this->ds->getType() == LNLDS)
  {
    VL(("Moreau::updateState -- LNLDS\n"));
    // get dynamical system
    LagrangianDS* d = static_cast<LagrangianDS*>(this->ds);
    // get velocity free, p, velocity and q pointers
    SimpleVector *vfree = d->getVelocityFreePtr();
    SimpleVector *p = d->getPPtr();
    SimpleVector *v = d->getVelocityPtr();
    SimpleVector *q = d->getQPtr();
    // temp temporary vector to save state k for velocity
    SimpleVector *temp = new SimpleVector((*v).size());
    *temp = *v ;
    // get time step
    double h = this->timeDiscretisation->getH();
    // compute state k+1 for velocity
    *v = *vfree + h * (*W) * (*p);
    // save state k in free state
    d->setQFree(d->getQ());
    d->setVelocityFree(*temp);
    delete temp;
    // get previous time step (i) state
    SimpleVector* vold, *qold;
    qold = static_cast<SimpleVector*>(d->getQMemoryPtr()->getSiconosVector(0));
    vold = static_cast<SimpleVector*>(d->getVelocityMemoryPtr()->getSiconosVector(0));
    // compute state k+1 for q
    *q = *qold + h * (theta * *v + (1.0 - theta)* *vold);
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
  this->W->display();
  cout << "| theta : " << this->theta << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}


void Moreau::saveIntegratorToXML()
{
  IN("Moreau::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (this->integratorxml != 0)
  {
    (static_cast<MoreauXML*>(this->integratorxml))->setTheta(this->theta);
    (static_cast<MoreauXML*>(this->integratorxml))->setW(this->W);
  }
  else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Moreau::saveIntegratorToXML\n");
}

void Moreau::saveWToXML()
{
  IN("Moreau::saveWToXML\n");
  if (this->integratorxml != 0)
  {
    (static_cast<MoreauXML*>(this->integratorxml))->setW(this->W);
  }
  else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Moreau::saveWToXML\n");
}



Moreau* Moreau::convert(OneStepIntegrator* osi)
{
  cout << "Moreau::convert (OneStepIntegrator* osi)" << endl;
  Moreau* moreau = dynamic_cast<Moreau*>(osi);
  return moreau;
}

// --- Default constructor ---
Moreau::Moreau(): OneStepIntegrator(), W(0), theta(0.1)
{
  this->integratorType = MOREAU_INTEGRATOR;
}
