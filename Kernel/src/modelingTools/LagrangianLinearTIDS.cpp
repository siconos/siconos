#include "LagrangianLinearTIDS.h"
using namespace std;

// --- Constructor from an xml file ---
LagrangianLinearTIDS::LagrangianLinearTIDS(DSXML * dsXML): LagrangianDS(dsXML), K(NULL), C(NULL),
  isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
  DSType = LTIDS;
  if (dsXML != NULL)
  {
    LagrangianLinearTIDSXML * lltidsxml = (static_cast <LagrangianLinearTIDSXML*>(dsxml));
    K = new SiconosMatrix(ndof, ndof);
    C = new SiconosMatrix(ndof, ndof);
    *K = lltidsxml->getK();
    *C = lltidsxml->getC();
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::LagrangianLinearTIDS - DSXML paramater must not be NULL");
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
}

// --- Constructor from a minimum set of data ---
LagrangianLinearTIDS::LagrangianLinearTIDS(const int& newNumber, const unsigned int& newNdof,
    const SimpleVector& newQ0, const SimpleVector& newVelocity0,
    const SiconosMatrix& newMass,
    const std::string& fExt, const SiconosMatrix& newK, const SiconosMatrix& newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, newMass), K(NULL), C(NULL), isKAllocatedIn(true), isCAllocatedIn(true)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  K = new SiconosMatrix(ndof, ndof);
  C = new SiconosMatrix(ndof, ndof);
  DSType = LTIDS;
  *K = newK;
  *C = newC;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}

LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
  if (isKAllocatedIn)
  {
    delete K;
    K = NULL;
  }
  if (isCAllocatedIn)
  {
    delete C;
    C = NULL;
  }
  OUT("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
}

void LagrangianLinearTIDS::setKPtr(SiconosMatrix *newPtr)
{
  if (isKAllocatedIn) delete K;
  K = newPtr;
  isKAllocatedIn = false;
}

void LagrangianLinearTIDS::setCPtr(SiconosMatrix *newPtr)
{
  if (isCAllocatedIn) delete C;
  C = newPtr;
  isCAllocatedIn = false;
}

void LagrangianLinearTIDS::display() const
{
  IN("LagrangianLinearTIDS::display\n");

  LagrangianDS::display();
  cout << "===== Lagrangian Linear Time Invariant System display ===== " << endl;
  cout << "- Stiffness Matrix K : " << endl;
  if (K != NULL) K->display();
  else cout << "-> NULL" << endl;
  cout << "- Viscosity Matrix C : " << endl;
  if (C != NULL) C->display();
  else cout << "-> NULL" << endl;
  cout << "=========================================================== " << endl;
  OUT("LagrangianLinearTIDS::display\n");
}

void LagrangianLinearTIDS::saveDSToXML()
{
  IN("LagrangianLinearTIDS::saveDSToXML\n");

  // --- general DS data---
  saveDSDataToXML();
  // --- other data  ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setNdof(ndof);
    lgptr->setMMatrix(mass);
    lgptr->setQ(q);
    lgptr->setQ0(q0);
    lgptr->setQMemory(qMemory);
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(fExtFunctionName);
    }
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setK(K);
    (static_cast <LagrangianLinearTIDSXML*>(dsxml))->setC(C);
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianLinearTIDS::saveDSToXML\n");
}

LagrangianLinearTIDS* LagrangianLinearTIDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianLinearTIDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianLinearTIDS* ltids = dynamic_cast<LagrangianLinearTIDS*>(ds);
  return ltids;
}

// --- Default (private) constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS(), K(NULL), C(NULL), isKAllocatedIn(false), isCAllocatedIn(false)
{
  DSType = LTIDS;
}

