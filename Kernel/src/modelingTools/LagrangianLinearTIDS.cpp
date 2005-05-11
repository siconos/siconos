#include "LagrangianLinearTIDS.h"
#include "check.h"


// --- Constructor from an xml file ---
LagrangianLinearTIDS::LagrangianLinearTIDS(DSXML * dsXML): LagrangianDS(dsXML), K(0), C(0)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
  this->DSType = LTIDS;
  if (dsXML != 0)
  {
    K = new SiconosMatrix(ndof, ndof);
    C = new SiconosMatrix(ndof, ndof);
    *this->K = (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->getK();
    *this->C = (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->getC();
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::LagrangianLinearTIDS - DSXML paramater must not be 0");
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
}

// --- Constructor from a minimum set of data ---
LagrangianLinearTIDS::LagrangianLinearTIDS(int newNumber, int newNdof,
    SiconosVector* newQ0, SiconosVector* newVelocity0, SiconosMatrix* newMass,
    string fExt, SiconosMatrix* newK, SiconosMatrix* newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, "BasicPlugin:computeMass", "BasicPlugin:computeFInt", fExt,
               "BasicPlugin:computeJacobianQFInt", "BasicPlugin:computeJacobianVelocityFInt",
               "BasicPlugin:computeJacobianQQNLInertia", "BasicPlugin:computeJacobianVelocityQNLInertia",
               "BasicPlugin:computeQNLInertia"), K(0), C(0)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  K = new SiconosMatrix(ndof, ndof);
  C = new SiconosMatrix(ndof, ndof);
  this->DSType = LTIDS;
  *mass = *newMass;
  *K = *newK;
  *C = *newC;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}


LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
  delete K;
  K = 0;
  delete C;
  C = 0;
  OUT("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
}

void LagrangianLinearTIDS::display() const
{
  IN("LagrangianLinearTIDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  LagrangianDS::display();
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LangrangianTIDS " << endl;
  cout << "| Stiffness Matrix K : " << endl;
  if (K != 0) K->display();
  else cout << "-> 0" << endl;
  cout << "| Viscosity Matrix C : " << endl;
  if (C != 0) C->display();
  else cout << "-> 0" << endl;
  cout << "-----------------------------------------------------" << endl << endl;
  OUT("LagrangianLinearTIDS::display\n");
}

void LagrangianLinearTIDS::saveDSToXML()
{
  IN("LagrangianLinearTIDS::saveDSToXML\n");

  // --- general DS data---
  saveDSDataToXML();
  // --- other data  ---
  if (this->dsxml != 0)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(this->dsxml);
    lgptr->setNdof(this->ndof);
    lgptr->setMMatrix(this->mass);
    lgptr->setQ(this->q);
    lgptr->setQ0(this->q0);
    lgptr->setQMemory(&(this->qMemory));
    lgptr->setVelocity(this->velocity);
    lgptr->setVelocity0(this->velocity0);
    lgptr->setVelocityMemory(&(this->velocityMemory));

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(this->fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(this->fExtFunctionName);
    }
    (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->setK(this->K);
    (static_cast <LagrangianLinearTIDSXML*>(this->dsxml))->setC(this->C);
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

// --- Default constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS(), K(0), C(0)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS() - Default constructor \n");
  this->DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS() - Default constructor \n");
}

