#include "LagrangianTIDS.h"

#include "check.h"


LagrangianTIDS::LagrangianTIDS()
{
  IN("LagrangianTIDS::LagrangianTIDS()\n");
  this->DSType = LTIDS;
  LagrangianNLDS::init();
  this->init();
  this->K = SiconosMatrix::SiconosMatrix();
  this->C = SiconosMatrix::SiconosMatrix();
  this->DSType = LTIDS;
  OUT("LagrangianTIDS::LagrangianTIDS()\n");
}

LagrangianTIDS::LagrangianTIDS(DSXML* dsxml): LagrangianNLDS(dsxml)
{
  IN("LagrangianTIDS::LagrangianTIDS(DSXML* dsxml):LagrangianNLDS(dsxml)\n");

  this->init();
  this->K = SiconosMatrix::SiconosMatrix();
  this->C = SiconosMatrix::SiconosMatrix();
  this->DSType = LTIDS;

  OUT("LagrangianTIDS::LagrangianTIDS(DSXML* dsxml):LagrangianNLDS(dsxml)\n");
}


LagrangianTIDS::~LagrangianTIDS()
{
  IN("LagrangianTIDS::~LagrangianTIDS()\n");

  OUT("LagrangianTIDS::~LagrangianTIDS()\n");
}

SiconosMatrix* LagrangianTIDS::getKPtr(void)
{
  return &this->K;
}

SiconosMatrix* LagrangianTIDS::getCPtr(void)
{
  return &this->C;
}


void LagrangianTIDS::init()
{
  IN("LagrangianTIDS::init\n");
  //LagrangianNLDS::init();
  this->computeFExtPtr = LTIDSComputeFExt;
  OUT("LagrangianTIDS::init\n");
}


void LagrangianTIDS::fillDSWithDSXML()
{
  IN("LagrangianTIDS::fillDSWithDSXML\n");

  if (this->dsxml != NULL)
  {
    this->K = (static_cast <LagrangianTIDSXML*>(this->dsxml))->getK();
    this->C = (static_cast <LagrangianTIDSXML*>(this->dsxml))->getC();
    //this->init();
    //    this->display();
  }
  else RuntimeException::selfThrow("LagrangianTIDS::fillDSWithDSXML - object DSXML does not exist");
  OUT("LagrangianTIDS::fillDSWithDSXML\n");
}

void LagrangianTIDS::display() const
{
  IN("LagrangianTIDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LangrangianTIDS " << endl;
  LagrangianNLDS::display();
  cout << "| Stiffness Matrix K : " << endl;
  K.display();
  cout << "| Viscosity Matrix C : " << endl;
  C.display();
  cout << "-----------------------------------------------------" << endl << endl;
  OUT("LagrangianTIDS::display\n");
}

void LagrangianTIDS::saveDSToXML()
{
  IN("LagrangianTIDS::saveDSToXML\n");
  LagrangianNLDS::saveDSToXML();

  if (this->dsxml != NULL)
  {
    (static_cast <LagrangianTIDSXML*>(this->dsxml))->setK(&(this->K));
    (static_cast <LagrangianTIDSXML*>(this->dsxml))->setC(&(this->C));
  }
  else RuntimeException::selfThrow("LagrangianTIDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianTIDS::saveDSToXML\n");
}


extern "C" double FextFunction(double time)
{
  double t = time;

  double res = -0.0;
  return res;
}

void LTIDSComputeFExt(int* sizeOfq, double* time, double* qPtr, double* fExt)
{
  IN("LagrangianTIDS : friend LTIDSComputeFInt\n");

  const double m = 1; // ball mass
  const double g = 9.8; // gravity

  cout << " /!\\!!! POUET !!!! [LagrangianTIDS : friend LTIDSComputeFInt]" << endl;
  getchar();

  int i;
  int n = *sizeOfq;
  double t = *time;

  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }

  fExt[0] = -m * g + FextFunction(t);
  OUT("LagrangianTIDS : friend LTIDSComputeFInt\n");
}

void LagrangianTIDS::createDynamicalSystem(DSXML * dsXML, int number, int ndof,
    SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
    string fExt, SiconosMatrix* K, SiconosMatrix* C) //,NSDS * nsds)
{
  if (dsXML != NULL)
  {
    LagrangianNLDS::createDynamicalSystem(dsXML);
    //this->init();
    LagrangianNLDS::fillDSWithDSXML();
    this->fillDSWithDSXML();
    this->linkDSXML();
    this->DSType = LTIDS;
  }
  else
  {
    this->number = number;
    this->ndof = ndof;

    LagrangianNLDS::SiconosVectorSizeInit();

    this->q0 = *q0;
    this->q = *q0;
    this->velocity0 = *velocity0;
    this->velocity = *velocity0;

    LagrangianNLDS::CompositeVectorInit();

    //this->setComputeMassFunction(this->cShared.getPluginName( mass ), this->cShared.getPluginFunctionName( mass ));
    this->mass = *mass;

    this->setComputeFExtFunction(this->cShared.getPluginName(fExt), this->cShared.getPluginFunctionName(fExt));

    //    this->setComputeJacobianQFIntFunction(this->cShared.getPluginName( jacobianQFInt ), this->cShared.getPluginFunctionName( jacobianQFInt ));
    //    this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName( jacobianVelocityFInt ), this->cShared.getPluginFunctionName( jacobianQFInt ));
    //    this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName( jacobianQQNLInertia ), this->cShared.getPluginFunctionName( jacobianQQNLInertia ));
    //    this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName( jacobianVelocityQNLInertia ), this->cShared.getPluginFunctionName( jacobianVelocityQNLInertia ));
    //
    //    this->setComputeQNLInertiaFunction(this->cShared.getPluginName( QNLInertia ), this->cShared.getPluginFunctionName( QNLInertia ));

    this->K = *K;
    this->C = *C;
    this->DSType = LTIDS;
  }
}


LagrangianTIDS* LagrangianTIDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianTIDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianTIDS* ltids = dynamic_cast<LagrangianTIDS*>(ds);
  return ltids;
}

