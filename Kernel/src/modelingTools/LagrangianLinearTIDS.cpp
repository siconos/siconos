#include "LagrangianLinearTIDS.h"
#include "check.h"


// --- Constructor from an xml file ---
LagrangianLinearTIDS::LagrangianLinearTIDS(DSXML * dsXML): LagrangianDS(dsXML), K(NULL), C(NULL)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
  DSType = LTIDS;
  if (dsXML != NULL)
  {
    K = new SiconosMatrix(ndof, ndof);
    C = new SiconosMatrix(ndof, ndof);
    *K = (static_cast <LagrangianLinearTIDSXML*>(dsxml))->getK();
    *C = (static_cast <LagrangianLinearTIDSXML*>(dsxml))->getC();
  }
  else RuntimeException::selfThrow("LagrangianLinearTIDS::LagrangianLinearTIDS - DSXML paramater must not be NULL");
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS() - Xml constructor\n");
}

// --- Constructor from a minimum set of data ---
LagrangianLinearTIDS::LagrangianLinearTIDS(int newNumber, int newNdof,
    SiconosVector* newQ0, SiconosVector* newVelocity0, SiconosMatrix* newMass,
    string fExt, SiconosMatrix* newK, SiconosMatrix* newC):
  LagrangianDS(newNumber, newNdof, newQ0, newVelocity0, "BasicPlugin:computeMass", "BasicPlugin:computeFInt", fExt,
               "BasicPlugin:computeJacobianQFInt", "BasicPlugin:computeJacobianVelocityFInt",
               "BasicPlugin:computeJacobianQQNLInertia", "BasicPlugin:computeJacobianVelocityQNLInertia",
               "BasicPlugin:computeQNLInertia"), K(NULL), C(NULL)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS -  Constructor from a minimum set of data\n");
  K = new SiconosMatrix(ndof, ndof);
  C = new SiconosMatrix(ndof, ndof);
  DSType = LTIDS;
  *mass = *newMass;
  *K = *newK;
  *C = *newC;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS - Constructor from a minimum set of data\n");
}


LagrangianLinearTIDS::~LagrangianLinearTIDS()
{
  IN("LagrangianLinearTIDS::~LagrangianLinearTIDS()\n");
  delete K;
  K = NULL;
  delete C;
  C = NULL;
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
  if (K != NULL) K->display();
  else cout << "-> NULL" << endl;
  cout << "| Viscosity Matrix C : " << endl;
  if (C != NULL) C->display();
  else cout << "-> NULL" << endl;
  cout << "-----------------------------------------------------" << endl << endl;
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
    lgptr->setQMemory(&(qMemory));
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(&(velocityMemory));

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

// --- Default constructor ---
LagrangianLinearTIDS::LagrangianLinearTIDS(): LagrangianDS(), K(NULL), C(NULL)
{
  IN("LagrangianLinearTIDS::LagrangianLinearTIDS() - Default constructor \n");
  DSType = LTIDS;
  OUT("LagrangianLinearTIDS::LagrangianLinearTIDS() - Default constructor \n");
}

