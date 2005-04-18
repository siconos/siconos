#include "LagrangianLinearTIDS.h"

#include "check.h"


LagrangianLinearTIDS::LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS()\n");
  this->DSType = LTIDS;
  LagrangianDS::init();
  this->init();
  this->K = SiconosMatrix::SiconosMatrix();
  this->C = SiconosMatrix::SiconosMatrix();
  this->DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS()\n");
}

LagrangianLinearTIDS::LagrangianLinearTIDS(DSXML * dsXML)//:LagrangianDS(dsXML)
{
  if (dsXML != NULL)
  {
    this->dsxml = dsXML;
    // /!\ not very good, to see if we can do better !!!
    LagrangianDS::LagrangianDS(dsXML);
    LagrangianDS::fillDSWithDSXML();
    this->fillDSWithDSXML();
    this->linkDSXML();
    this->DSType = LTIDS;
  }
  else
  {
    cout << "LagrangianLinearTIDS::LagrangianLinearTIDS - DSXML paramater musn't be NULL" << endl;
  }
}

LagrangianLinearTIDS::LagrangianLinearTIDS(int number, int ndof,
    SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
    string fExt, SiconosMatrix* K, SiconosMatrix* C)
{
  this->number = number;
  this->ndof = ndof;

  LagrangianDS::SiconosVectorSizeInit();

  this->q0 = *q0;
  this->q = *q0;
  this->velocity0 = *velocity0;
  this->velocity = *velocity0;

  LagrangianDS::CompositeVectorInit();

  this->mass = *mass;

  this->setComputeFExtFunction(this->cShared.getPluginName(fExt), this->cShared.getPluginFunctionName(fExt));
  this->K = *K;
  this->C = *C;
  this->DSType = LTIDS;
}


LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");

  OUT("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
}

SiconosMatrix* LagrangianLinearTIDS::getKPtr(void)
{
  return &this->K;
}

SiconosMatrix* LagrangianLinearTIDS::getCPtr(void)
{
  return &this->C;
}


void LagrangianLinearTIDS::init()
{
  IN("LagrangianLinearTIDS::init\n");
  //LagrangianDS::init();
  this->computeFExtPtr = LTIDSComputeFExt;
  OUT("LagrangianLinearTIDS::init\n");
}


void LagrangianLinearTIDS::fillDSWithDSXML()
{
  IN("LagrangianLinearTIDS::fillDSWithDSXML\n");

  if (this->dsxml != NULL)
  {
    this->K = (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->getK();
    this->C = (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->getC();
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::fillDSWithDSXML - object DSXML does not exist");
  OUT("LagrangianLinearTIDS::fillDSWithDSXML\n");
}

void LagrangianLinearTIDS::display() const
{
  IN("LagrangianLinearTIDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LangrangianTIDS " << endl;
  LagrangianDS::display();
  cout << "| Stiffness Matrix K : " << endl;
  K.display();
  cout << "| Viscosity Matrix C : " << endl;
  C.display();
  cout << "-----------------------------------------------------" << endl << endl;
  OUT("LagrangianLinearTIDS::display\n");
}

void LagrangianLinearTIDS::saveDSToXML()
{
  IN("LagrangianLinearTIDS::saveDSToXML\n");
  LagrangianDS::saveDSToXML();

  if (this->dsxml != NULL)
  {
    (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->setK(&(this->K));
    (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->setC(&(this->C));
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianLinearTIDS::saveDSToXML\n");
}


extern "C" double FextFunction(double time)
{
  double t = time;

  double res = -0.0;
  return res;
}

void LTIDSComputeFExt(int* sizeOfq, double* time, double* qPtr, double* fExt)
{
  IN("LagrangianLinearTIDS : friend LTIDSComputeFInt\n");

  const double m = 1; // ball mass
  const double g = 9.8; // gravity

  cout << " /!\\!!! POUET !!!! [LagrangianLinearTIDS : friend LTIDSComputeFInt]" << endl;
  getchar();

  int i;
  int n = *sizeOfq;
  double t = *time;

  for (i = 0; i < n; i++)
  {
    fExt[i] = 0.0;
  }

  fExt[0] = -m * g + FextFunction(t);
  OUT("LagrangianLinearTIDS : friend LTIDSComputeFInt\n");
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianLinearTIDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

