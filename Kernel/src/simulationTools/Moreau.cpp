#include "Moreau.h"
// includes to be deleted thanks to factories:
#include "LagrangianLinearTIDS.h"
#include "LagrangianDS.h"
#include "LinearDS.h"

using namespace std;


// --- xml constructor ---
Moreau::Moreau(OneStepIntegratorXML *osiXML, TimeDiscretisation* td, DynamicalSystem* ds) :
  OneStepIntegrator(td, ds), W(NULL), isWAllocatedIn(true), theta(0.1)
{

  integratorXml = osiXML;
  integratorType = MOREAU_INTEGRATOR;
  // Memory allocation for W
  int sizeW = (static_cast<LagrangianDS*>(ds))->getQPtr()->size();
  W = new SiconosMatrix(sizeW, sizeW);

  // xml loading
  if (osiXML != NULL)
  {
    if (integratorXml->hasR()) sizeMem = integratorXml->getR();
    if ((static_cast<MoreauXML*>(integratorXml))->hasW() == true)
    {
      *W = (static_cast<MoreauXML*>(integratorXml))->getW();
    }
    if ((static_cast<MoreauXML*>(integratorXml))->hasTheta() == true)
    {
      theta = (static_cast<MoreauXML*>(integratorXml))->getTheta();
    }
  }
  else RuntimeException::selfThrow("Moreau::Moreau() - xml constructor - IntegratorXML object not exists");
}

// --- constructor from a minimum set of data ---
Moreau::Moreau(TimeDiscretisation* td, DynamicalSystem* ds, const double& newTheta):
  OneStepIntegrator(td, ds), W(NULL), isWAllocatedIn(false), theta(newTheta)
{
  integratorType = MOREAU_INTEGRATOR;
  // Memory allocation for W
  string type = ds->getType();
  if (type == LNLDS || type == LTIDS)
  {
    int sizeW = (static_cast<LagrangianDS*>(ds))->getQPtr()->size();
    W = new SiconosMatrix(sizeW, sizeW);
    isWAllocatedIn = true;
  }
  else if (type == LDS)
  {
    int sizeW = (static_cast<LinearDS*>(ds))->getN();
    W = new SiconosMatrix(sizeW, sizeW);
    isWAllocatedIn = true;
  }
  else
    RuntimeException::selfThrow("Moreau::Moreau() - constructor from data - Not yet implemented for ds type " + type);
}

Moreau::~Moreau()
{
  if (isWAllocatedIn)
  {
    delete W;
    W = NULL;
  }
  if (ds != NULL && ds->getType() == LNLDS) ds->freeTmpWorkVector("LagNLDSMoreau");
}

void Moreau::initialize()
{
  IN("Moreau::initialize\n");
  OneStepIntegrator::initialize();
  // Get initial time
  double t0 = timeDiscretisation->getT0();
  // Compute W
  computeW(t0);
  if (ds->getType() == LNLDS)
    ds->allocateTmpWorkVector("LagNLDSMoreau", W->size(0));
  OUT("Moreau::initialize\n");
}

void Moreau::computeW(const double& t)
{
  double h = timeDiscretisation->getH(); // time step

  // === Lagrangian dynamical system
  if (ds->getType() == LNLDS)
  {
    LagrangianDS* d = static_cast<LagrangianDS*>(ds);
    // Check if W is allocated; if not, do allocation.
    if (W == NULL)
    {
      unsigned int size = d->getQPtr()->size();
      W = new SiconosMatrix(size, size);
      isWAllocatedIn = true;
    }
    // Compute Mass matrix (if loaded from plugin)
    if (d->getIsLDSPlugin(0))
      d->computeMass(t);
    // Compute and get Jacobian (if loaded from plugin)
    if (d->getIsLDSPlugin(4))
      d->computeJacobianQFInt(t);
    if (d->getIsLDSPlugin(5))
      d->computeJacobianVelocityFInt(t);
    if (d->getIsLDSPlugin(6))
      d->computeJacobianQNNL();
    if (d->getIsLDSPlugin(7))
      d->computeJacobianVelocityNNL();

    SiconosMatrix *KFint, *KQNL, *CFint, *CQNL ;
    KFint = d->getJacobianQFIntPtr();
    KQNL  = d->getJacobianQNNLPtr();
    CFint = d->getJacobianVelocityFIntPtr();
    CQNL  = d->getJacobianVelocityNNLPtr();

    // Get Mass matrix
    SiconosMatrix *M = d->getMassPtr();

    // Compute W
    *W = *M + h * theta * (*CFint + *CQNL  + h * theta * (*KFint + *KQNL));
  }
  // === Lagrangian linear time invariant system ===
  else if (ds->getType() == LTIDS)
  {
    LagrangianDS* d = static_cast<LagrangianDS*>(ds);
    // Check if W is allocated; if not, do allocation.
    if (W == NULL)
    {
      unsigned int size = d->getQPtr()->size();
      W = new SiconosMatrix(size, size);
      isWAllocatedIn = true;
    }
    // Get K, C and Mass
    SiconosMatrix *M, *K, *C ;
    K = ((static_cast<LagrangianLinearTIDS*>(d))->getKPtr());
    C = ((static_cast<LagrangianLinearTIDS*>(d))->getCPtr());
    M = d->getMassPtr();

    // Compute W
    *W = *M + h * theta * (*C + h * theta* *K);
  }

  // === Linear dynamical system ===
  else if (ds->getType() == LDS)
  {
    LinearDS* d = static_cast<LinearDS*>(ds);
    SiconosMatrix *I;
    int size = d->getN();
    // Check if W is allocated
    if (W == NULL)
    {
      W = new SiconosMatrix(size, size);
      isWAllocatedIn = true;
    }
    I = new SiconosMatrix(size, size);
    I->eye();
    *W = *I - (h * theta * (d->getA()));
    delete I;
  }
  // === ===
  else RuntimeException::selfThrow("Moreau::computeW - not yet implemented for Dynamical system type :" + ds->getType());
  // LU factorization of W
  W->PLUFactorizationInPlace();
  // At the time, W inverse is saved in Moreau object -> \todo: to be reviewed: use forwarBackward in OneStepNS formalize to avoid inversion of W => work on Mlcp
  W->PLUInverseInPlace();
}


void Moreau::computeFreeState()
{
  IN("Moreau::computeFreeState\n");
  // get current time, theta and time step
  double t = timeDiscretisation->getStrategyPtr()->getModelPtr()->getCurrentT();
  double h = timeDiscretisation->getH();
  // Previous time step (i)
  double told = t - h;

  // Get the DS type
  std::string dstyp = ds->getType();

  if ((dstyp == LNLDS) || (dstyp == LTIDS))
  {

    // -- Get the DS --
    LagrangianDS* d = static_cast<LagrangianDS*>(ds);

    // --- RESfree calculus ---
    //
    // Get state i (previous time step)
    SimpleVector* vold, *qold;
    qold = static_cast<SimpleVector*>(d->getQMemoryPtr()->getSiconosVector(0));
    vold = static_cast<SimpleVector*>(d->getVelocityMemoryPtr()->getSiconosVector(0));

    // Computation of the external forces
    d->computeFExt(told);
    SimpleVector FExt0 = d->getFExt();
    d->computeFExt(t);
    SimpleVector FExt1 = d->getFExt();
    // RESfree ...
    SimpleVector *v = d->getVelocityPtr();
    SimpleVector *RESfree = new SimpleVector(FExt1.size());
    // Velocity free
    SimpleVector *vfree = d->getVelocityFreePtr();

    // --- Compute Velocity Free ---
    // For general Lagrangian system LNLDS:
    if (dstyp == LNLDS)
    {
      // Get Mass (remark: M is computed for present state during computeW(t) )
      SiconosMatrix *M = d -> getMassPtr();
      // Compute Qint and Fint
      // for state i
      // warning: get values and not pointers
      d->computeNNL(qold, vold);
      d->computeFInt(told, qold, vold);
      SimpleVector QNL0 = d->getNNL();
      SimpleVector FInt0 = d->getFInt();
      // for present state
      // warning: get values and not pointers
      d->computeNNL();
      d->computeFInt(t);
      SimpleVector QNL1 = d->getNNL();
      SimpleVector FInt1 = d->getFInt();
      // Compute ResFree and vfree solution of Wk(v-vfree)=RESfree
      *RESfree = *M * (*v - *vold) + h * ((1.0 - theta) * (QNL0 + FInt0 - FExt0) + theta * (QNL1 + FInt1 - FExt1));
      *vfree = *v - *W * *RESfree;
    }
    // --- For linear Lagrangian LTIDS:
    else
    {
      // get K, M and C mass pointers
      SiconosMatrix * K = static_cast<LagrangianLinearTIDS*>(d)->getKPtr();
      SiconosMatrix * C = static_cast<LagrangianLinearTIDS*>(d)->getCPtr();
      // Compute ResFree and vfree
      *RESfree = -h * (theta * FExt1 + (1.0 - theta) * FExt0 - (*C * *vold) - (*K * *qold) - h * theta * (*K * *vold));
      *vfree =  *vold - *W * *RESfree;
    }
    // calculate qfree (whereas it is useless for future computation)
    SimpleVector *qfree = d->getQFreePtr();
    *qfree = (*qold) + h * (theta * (*vfree) + (1.0 - theta) * (*vold));
    delete RESfree;
  }
  else if (dstyp == LDS)
  {
    SiconosVector *xfree = ds->getXFreePtr();
    SiconosVector *xold = ds->getXMemoryPtr()->getSiconosVector(0);
    if (theta == 1.0)
      *xfree = *W * (*xold);
    else
    {
      unsigned int size = ds->getN();
      SiconosMatrix *I = new SiconosMatrix(size, size);
      I->eye();
      SiconosMatrix * A = static_cast<LinearDS*>(ds)->getAPtr();
      *xfree = *W * (*I + h* *A * (1.0 - theta))* *xold;
      delete I;
    }
  }
  else RuntimeException::selfThrow("Moreau::computeFreeState - not yet implemented for Dynamical system type: " + dstyp);
  OUT("Moreau::computeFreeState\n");
}


void Moreau::integrate()
{
  IN("Moreau::integrate()\n");

  double h = timeDiscretisation->getH();
  double t = timeDiscretisation->getStrategyPtr()->getModelPtr()->getCurrentT();
  double told = t - h;

  if (ds->getType() == LNLDS)
  {
    //VL(("Moreau::integrate -- LNLDS\n"));
    // We do not use integrate() for LNDS
  }
  else if (ds->getType() == LTIDS)
  {
    VL(("Moreau::integrate -- LTIDS\n"));
    // get the ds
    LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(ds);
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
    *v = *vold + *W * (h * (theta * FExt1 + (1.0 - theta) * FExt0 - (*C * *vold) - (*K * *qold) - h * theta * (*K * *vold)) + *p);
    // q computation
    *q = (*qold) + h * ((theta * (*v)) + (1.0 - theta) * (*vold));
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

  double h = timeDiscretisation->getH();

  // Get the DS type
  std::string dsType = ds->getType();

  if ((dsType == LNLDS) || (dsType == LTIDS))
  {
    // get dynamical system
    LagrangianDS* d = static_cast<LagrangianDS*>(ds);
    // get velocity free, p, velocity and q pointers
    SimpleVector *vfree = d->getVelocityFreePtr();
    SimpleVector *p = d->getPPtr();
    SimpleVector *v = d->getVelocityPtr();
    SimpleVector *q = d->getQPtr();
    // Save value of q and v in stateTmp for future convergence computation
    if (dsType == LNLDS)
      ds->addTmpWorkVector(v, "LagNLDSMoreau");
    // Compute velocity
    *v = *vfree +  *W * *p;
    // Compute q
    //  -> get previous time step state
    SimpleVector *vold = static_cast<SimpleVector*>(d->getVelocityMemoryPtr()->getSiconosVector(0));
    SimpleVector *qold = static_cast<SimpleVector*>(d->getQMemoryPtr()->getSiconosVector(0));
    *q = *qold + h * (theta * *v + (1.0 - theta)* *vold);
    // set reaction to zero
    p->zero();
    // --- Update W for general Lagrangian system
    if (dsType == LNLDS)
    {
      double t = timeDiscretisation->getStrategyPtr()->getModelPtr()->getCurrentT();
      computeW(t);
    }
    // Remark: for Linear system, W is already saved in object member w
  }
  else if (dsType == LDS)
  {
    SiconosVector* x = ds->getXPtr();
    SiconosVector* xFree = ds->getXFreePtr();

    *x = *xFree + h * *W * *(ds->getRPtr()) ;
  }
  else RuntimeException::selfThrow("Moreau::updateState - not yet implemented for Dynamical system type: " + dsType);
  OUT("Moreau::updateState\n");
}


void Moreau::display() const
{
  OneStepIntegrator::display();

  cout << "====== Moreau OSI display ======" << endl;
  cout << "____ data of the Moreau Integrator " << endl;
  cout << "| W " << endl;
  if (W != NULL) W->display();
  else cout << "-> NULL" << endl;
  cout << "| theta : " << theta << endl;
  cout << "================================" << endl;
}

void Moreau::saveIntegratorToXML()
{
  IN("Moreau::saveIntegratorToXML\n");
  OneStepIntegrator::saveIntegratorToXML();
  if (integratorXml != NULL)
  {
    (static_cast<MoreauXML*>(integratorXml))->setTheta(theta);
    (static_cast<MoreauXML*>(integratorXml))->setW(W);
  }
  else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  OUT("Moreau::saveIntegratorToXML\n");
}

void Moreau::saveWToXML()
{
  IN("Moreau::saveWToXML\n");
  if (integratorXml != NULL)
  {
    (static_cast<MoreauXML*>(integratorXml))->setW(W);
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
Moreau::Moreau(): OneStepIntegrator(), W(NULL), isWAllocatedIn(false), theta(0.1)
{
  integratorType = MOREAU_INTEGRATOR;
}
