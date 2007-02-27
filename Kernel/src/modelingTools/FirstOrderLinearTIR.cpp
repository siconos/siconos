/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "FirstOrderLinearTIR.h"
#include "FirstOrderLinearTIRXML.h"
#include "Interaction.h"

using namespace std;

// Default (private) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(): FirstOrderR("LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  isAllocatedIn["C"] = false;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
  isAllocatedIn["e"] = false;
  isAllocatedIn["B"] = false;
}

// xml constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(RelationXML* relxml):
  FirstOrderR(relxml, "LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  FirstOrderLinearTIRXML * lTIRxml = static_cast<FirstOrderLinearTIRXML *>(relationxml);
  // get matrices values.
  // The only required variables are C and B, other are optional.
  // If ouput plug-in is provided ...
  if (isPlugged["output"] || isPlugged["input"])
    RuntimeException::selfThrow("FirstOrderLinearTIR xml constructor failed, you can not give plug-in functions for input and output.");

  if (!lTIRxml->hasC())
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor: input matrix C is missing in xml file ");

  C = new SimpleMatrix(lTIRxml->getC());
  isAllocatedIn["C"] = true;

  if (!lTIRxml->hasB())
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor: input matrix B is missing in xml file ");

  B = new SimpleMatrix(lTIRxml->getB());
  isAllocatedIn["B"] = true;

  isAllocatedIn["D"] = false;
  if (lTIRxml->hasD())
  {
    D = new SimpleMatrix(lTIRxml->getD());
    isAllocatedIn["D"] = true;
  }

  isAllocatedIn["F"] = false;
  if (lTIRxml->hasF())
  {
    F = new SimpleMatrix(lTIRxml->getF());
    isAllocatedIn["F"] = true;
  }

  isAllocatedIn["e"] = false;
  if (lTIRxml->hasE())
  {
    e = new SimpleVector(lTIRxml->getE());
    isAllocatedIn["e"] = true;
  }
}

// Minimum data (C, B) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderR("LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  C = new SimpleMatrix(newC);
  B = new SimpleMatrix(newB);
  isAllocatedIn["C"] = true;
  isAllocatedIn["B"] = true;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
  isAllocatedIn["e"] = false;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
    const SiconosMatrix& newF, const SimpleVector& newE,
    const SiconosMatrix& newB):
  FirstOrderR("LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  C = new SimpleMatrix(newC);
  isAllocatedIn["C"] = true;

  D = new SimpleMatrix(newD);
  isAllocatedIn["D"] = true;

  F = new SimpleMatrix(newF);
  isAllocatedIn["F"] = true;

  e = new SimpleVector(newE);
  isAllocatedIn["e"] = true;

  B = new SimpleMatrix(newB);
  isAllocatedIn["B"] = true;
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SiconosMatrix * newC, SiconosMatrix * newB):
  FirstOrderR("LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  isAllocatedIn["C"] = false;
  isAllocatedIn["B"] = false;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
  isAllocatedIn["e"] = false;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SiconosMatrix* newC, SiconosMatrix* newD, SiconosMatrix* newF, SimpleVector* newE, SiconosMatrix* newB):
  FirstOrderR("LinearTIR"), C(NULL), D(NULL), F(NULL), e(NULL), B(NULL)
{
  isAllocatedIn["C"] = false;
  isAllocatedIn["B"] = false;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
  isAllocatedIn["e"] = false;
}

FirstOrderLinearTIR::~FirstOrderLinearTIR()
{
  if (isAllocatedIn["C"]) delete C;
  C = NULL;
  if (isAllocatedIn["D"]) delete D;
  D = NULL;
  if (isAllocatedIn["F"]) delete F;
  F = NULL;
  if (isAllocatedIn["e"]) delete e;
  e = NULL;
  if (isAllocatedIn["B"]) delete B;
  B = NULL;
}

void FirstOrderLinearTIR::initialize()
{
  FirstOrderR::initialize();

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = interaction->getInteractionSize();
  unsigned int sizeX = interaction->getSizeOfDS();

  if (C->size(0) != sizeY || C->size(1) != sizeX)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  if (B->size(0) != sizeX || B->size(1) != sizeY)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize , inconsistent size between C and B.");

  if (D != NULL)
  {
    if (D->size(0) != sizeY || D->size(1) != sizeY)
      RuntimeException::selfThrow("FirstOrderLinearTIR::initialize , inconsistent size between C and D.");
  }

  if (F != NULL && (F->size(0) != sizeY))
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize , inconsistent size between C and F.");

  if (e != NULL && e->size() != sizeY)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize , inconsistent size between C and e.");
}

// setters

void FirstOrderLinearTIR::setC(const SiconosMatrix& newValue)
{
  if (C == NULL)
  {
    C = new SimpleMatrix(newValue);
    isAllocatedIn["C"] = true;
  }
  else
  {
    if ((newValue.size(0) == C->size(0)) && (newValue.size(1) == C->size(1)))
      *C = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearTIR - setC: inconsistent dimensions with problem size for input matrix C.");
  }
}

void FirstOrderLinearTIR::setCPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["C"]) delete C;
  C = newPtr;
  isAllocatedIn["C"] = false;
}

void FirstOrderLinearTIR::setD(const SiconosMatrix& newValue)
{
  if (D == NULL)
  {
    D = new SimpleMatrix(newValue);
    isAllocatedIn["D"] = true;
  }
  else
  {
    if ((newValue.size(0) == D->size(0)) && (newValue.size(1) == D->size(1)))
      *D = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearTIR - setD: inconsistent dimensions with problem size for input matrix D.");
  }
}

void FirstOrderLinearTIR::setDPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["D"]) delete D;
  D = newPtr;
  isAllocatedIn["D"] = false;
}

void FirstOrderLinearTIR::setF(const SiconosMatrix& newValue)
{
  if (F == NULL)
  {
    F = new SimpleMatrix(newValue);
    isAllocatedIn["F"] = true;
  }
  else
  {
    if ((newValue.size(0) == F->size(0)) && (newValue.size(1) == F->size(1)))
      *F = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearTIR - setF: inconsistent dimensions with problem size for input matrix F.");
  }
}

void FirstOrderLinearTIR::setFPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["F"]) delete F;
  F = newPtr;
  isAllocatedIn["F"] = false;
}

void FirstOrderLinearTIR::setE(const SiconosVector& newValue)
{
  if (e == NULL)
  {
    e = new SimpleVector(newValue);
    isAllocatedIn["e"] = true;
  }
  else
  {
    if (newValue.size() == e->size())
      *e = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearTIR - setE: inconsistent dimensions with problem size for input vector e.");
  }
}

void FirstOrderLinearTIR::setEPtr(SiconosVector* newPtr)
{
  if (isAllocatedIn["e"]) delete e;
  e = newPtr;
  isAllocatedIn["e"] = false;
}

void FirstOrderLinearTIR::setB(const SiconosMatrix& newValue)
{
  if (B == NULL)
  {
    B = new SimpleMatrix(newValue);
    isAllocatedIn["B"] = true;
  }
  else
  {
    if ((newValue.size(0) == B->size(0)) && (newValue.size(1) == B->size(1)))
      *B = newValue;
    else
      RuntimeException::selfThrow("FirstOrderLinearTIR - setB: inconsistent dimensions with problem size for input matrix B.");
  }
}

void FirstOrderLinearTIR::setBPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["B"]) delete B;
  B = newPtr;
  isAllocatedIn["B"] = false;
}

void FirstOrderLinearTIR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SiconosVector *y = interaction->getYPtr(0);
  SiconosVector *lambda = interaction->getLambdaPtr(0);

  // compute y
  *y = prod(*C, *data["x"]);

  if (D != NULL)
    *y += prod(*D, *lambda);

  if (e != NULL)
    *y += *e;

  if (F != NULL)
    *y += prod(*F, *data["z"]);
}

void FirstOrderLinearTIR::computeFreeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  SiconosVector *yFree = interaction->getYPtr(0);
  // warning : yFree is saved in y !!

  // compute yFree
  *yFree = prod(*C, *data["xFree"]);

  if (e != NULL)
    *yFree += *e;

  if (F != NULL)
    *yFree += prod(*F, *data["z"]);
}

void FirstOrderLinearTIR::computeInput(double time, unsigned int level)
{
  // We get lambda of the interaction (pointers)
  SiconosVector *lambda = interaction->getLambdaPtr(level);
  *data["r"] += prod(*B, *lambda);
}

void FirstOrderLinearTIR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (C != NULL) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D != NULL) D->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F != NULL) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e != NULL) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (B != NULL) B->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearTIR::saveRelationToXML() const
{
  if (relationxml == NULL)
    RuntimeException::selfThrow("FirstOrderLinearTIR::saveRelationToXML, no xml object found.");

  FirstOrderLinearTIRXML * lTIRxml = (static_cast<FirstOrderLinearTIRXML*>(relationxml));
  lTIRxml->setC(*C);
  lTIRxml->setD(*D);
  lTIRxml->setF(*F);
  lTIRxml->setE(*e);
  lTIRxml->setB(*B);
}

FirstOrderLinearTIR* FirstOrderLinearTIR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearTIR*>(r);
}

