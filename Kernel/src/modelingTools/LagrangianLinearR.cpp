/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearR.h"
#include "LagrangianLinearRXML.h"
#include "Interaction.h"
//
#include "LagrangianDS.h"

using namespace std;

void LagrangianLinearR::initAllocationFlags(bool in)
{
  isAllocatedIn["H"] = in;
  isAllocatedIn["b"] = in;
  isAllocatedIn["D"] = in;
  isAllocatedIn["F"] = in;
}

// Default (private) constructor
LagrangianLinearR::LagrangianLinearR():
  LagrangianR("LinearR"), H(NULL), b(NULL), D(NULL), F(NULL)
{
  initAllocationFlags(false);
}

// Xml constructor
LagrangianLinearR::LagrangianLinearR(RelationXML* relxml):
  LagrangianR(relxml, "LinearR"), H(NULL), b(NULL), D(NULL), F(NULL)
{
  LagrangianLinearRXML* LLRxml = (static_cast<LagrangianLinearRXML*>(relationxml));

  initAllocationFlags(false);

  // H is the minimal required input.
  if (!LLRxml->hasH())
    RuntimeException::selfThrow("LagrangianLinearR:: xml constructor failed: can not find input for H matrix.");

  H = new SimpleMatrix(LLRxml->getH());
  isAllocatedIn["H"] = true;

  if (LLRxml->hasB())
  {
    b = new SimpleVector(LLRxml->getB());
    isAllocatedIn["b"] = true;
  }

  if (LLRxml->hasD())
  {
    D = new SimpleMatrix(LLRxml->getD());
    isAllocatedIn["D"] = true;
  }

  if (LLRxml->hasF())
  {
    F = new SimpleMatrix(LLRxml->getF());
    isAllocatedIn["F"] = true;
  }
}

// Constructor from data: H, b.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB):
  LagrangianR("LinearR"), H(NULL), b(NULL), D(NULL), F(NULL)
{
  H = new SimpleMatrix(newH);
  isAllocatedIn["H"] = true;

  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
}

// Constructor from data: H.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH):
  LagrangianR("LinearR"), H(NULL), b(NULL), D(NULL), F(NULL)
{
  H = new SimpleMatrix(newH);
  isAllocatedIn["H"] = true;
  isAllocatedIn["b"] = false;
  isAllocatedIn["D"] = false;
  isAllocatedIn["F"] = false;
}

// Constructor from data: H, b, and D.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB, const SiconosMatrix& newD):
  LagrangianR("LinearR"), H(NULL), b(NULL), D(NULL)
{
  H = new SimpleMatrix(newH);
  isAllocatedIn["H"] = true;
  D = new SimpleMatrix(newD);
  isAllocatedIn["D"] = true;
  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  isAllocatedIn["F"] = false;
}

// Constructor from data: H, b, D and F.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB, const SiconosMatrix& newD, const SiconosMatrix& newF):
  LagrangianR("LinearR"), H(NULL), b(NULL), D(NULL)
{
  H = new SimpleMatrix(newH);
  isAllocatedIn["H"] = true;
  D = new SimpleMatrix(newD);
  isAllocatedIn["D"] = true;
  b = new SimpleVector(newB);
  isAllocatedIn["b"] = true;
  F = new SimpleMatrix(newF);
  isAllocatedIn["F"] = true;
}

LagrangianLinearR::~LagrangianLinearR()
{
  if (isAllocatedIn["H"]) delete H;
  H = NULL;
  if (isAllocatedIn["b"]) delete b;
  b = NULL;
  if (isAllocatedIn["D"]) delete D;
  D = NULL;
  if (isAllocatedIn["F"]) delete F;
  F = NULL;
}

void LagrangianLinearR::initComponents()
{
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeDS = interaction->getSizeOfDS();

  if (H->size(1) != sizeDS || H->size(0) != sizeY)
    RuntimeException::selfThrow("LagrangianLinearR::initComponents inconsistent sizes between H matrix and the interaction.");

  if (D != NULL)
    if (D->size(0) != sizeY || D->size(1) != sizeY)
      RuntimeException::selfThrow("LagrangianLinearR::initComponents inconsistent sizes between D matrix and the interaction.");
  if (b != NULL)
    if (b->size() != sizeY)
      RuntimeException::selfThrow("LagrangianLinearR::initComponents inconsistent sizes between b vector and the dimension of the interaction.");

  if (F != NULL)
  {
    unsigned int sizeZ = interaction->getSizeZ();
    if (F->size(0) != sizeY || F->size(1) != sizeZ)
      RuntimeException::selfThrow("LagrangianLinearR::initComponents inconsistent sizes between F matrix and the interaction.");
  }

}

// Setters

void LagrangianLinearR::setH(const SiconosMatrix& newValue)
{
  if (H == NULL)
  {
    H = new SimpleMatrix(newValue);
    isAllocatedIn["H"] = true;
  }
  else
    *H = newValue;
}

void LagrangianLinearR::setHPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["H"]) delete H;
  H = newPtr;
  isAllocatedIn["H"] = false;
}

void LagrangianLinearR::setB(const SimpleVector& newValue)
{
  if (b == NULL)
  {
    b = new SimpleVector(newValue);
    isAllocatedIn["b"] = true;
  }
  else
    *b = newValue;
}

void LagrangianLinearR::setBPtr(SimpleVector *newPtr)
{
  if (isAllocatedIn["b"]) delete b;
  b = newPtr;
  isAllocatedIn["b"] = false;
}

void LagrangianLinearR::setD(const SiconosMatrix& newValue)
{
  if (D == NULL)
  {
    D = new SimpleMatrix(newValue);
    isAllocatedIn["D"] = true;
  }
  else
    *D = newValue;
}

void LagrangianLinearR::setDPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["D"]) delete D;
  D = newPtr;
  isAllocatedIn["D"] = false;
}

void LagrangianLinearR::setF(const SiconosMatrix& newValue)
{
  if (F == NULL)
  {
    F = new SimpleMatrix(newValue);
    isAllocatedIn["F"] = true;
  }
  else
    *F = newValue;
}

void LagrangianLinearR::setFPtr(SiconosMatrix *newPtr)
{
  if (isAllocatedIn["F"]) delete F;
  F = newPtr;
  isAllocatedIn["F"] = false;
}

void LagrangianLinearR::computeOutput(double time, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SiconosVector *y = interaction->getYPtr(derivativeNumber);
  SiconosVector *lambda = interaction->getLambdaPtr(derivativeNumber);

  string name = "q" + toString<unsigned int>(derivativeNumber);

  prod(*H, *data[name], *y);
  if (derivativeNumber == 0 && b != NULL)
    *y += *b;

  if (D != NULL)
    *y += prod(*D, *lambda) ;

  if (F != NULL)
    *y += prod(*F, *data["z"]);
}

void LagrangianLinearR::computeFreeOutput(double time, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SiconosVector *y = interaction->getYPtr(derivativeNumber);
  SiconosVector *lambda = interaction->getLambdaPtr(derivativeNumber);

  string name = "q" + toString<unsigned int>(derivativeNumber) + "Free";

  if (derivativeNumber == 2) name = "q2";

  prod(*H, *data[name], *y);
  if (derivativeNumber == 0 && b != NULL)
    *y += *b;

  if (D != NULL)
    *y += prod(*D, *lambda) ;

  if (F != NULL)
    *y += prod(*F, *data["z"]);
}

void LagrangianLinearR::computeInput(double time, const unsigned int level)
{
  // get lambda of the concerned interaction
  string name = "p" + toString<unsigned int>(level);

  SiconosVector *lambda = new SimpleVector(*interaction->getLambdaPtr(level));
  // compute p = Ht lambda
  SiconosMatrix * HT = new SimpleMatrix(*H);

  HT->trans();
  *data[name] += prod(*HT, *lambda);

  delete HT;
  delete lambda;
  //gemv(CblasTrans,1.0,*H,*lambda,1.0, *data[name]); => not yet implemented for BlockVectors.
}

void LagrangianLinearR::saveRelationToXML() const
{
  if (relationxml == NULL)
    RuntimeException::selfThrow("LagrangianLinearR::saveRelationToXML - object RelationXML does not exist");

  (static_cast<LagrangianLinearRXML*>(relationxml))->setH(*H) ;
  (static_cast<LagrangianLinearRXML*>(relationxml))->setB(*b) ;
  (static_cast<LagrangianLinearRXML*>(relationxml))->setD(*D) ;
  (static_cast<LagrangianLinearRXML*>(relationxml))->setF(*F) ;
}

LagrangianLinearR* LagrangianLinearR::convert(Relation *r)
{
  return dynamic_cast<LagrangianLinearR*>(r);
}

void LagrangianLinearR::display() const
{
  LagrangianR::display();
  cout << "===== Lagrangian Linear Relation display ===== " << endl;
  cout << " H: " << endl;
  if (H != NULL)
    H->display();
  else
    cout << " -> NULL " << endl;
  cout << " b: " << endl;
  if (b != NULL)
    b->display();
  else
    cout << " -> NULL " << endl;
  cout << " D: " << endl;
  if (D != NULL)
    D->display();
  else
    cout << " -> NULL " << endl;
  cout << " F: " << endl;
  if (F != NULL)
    F->display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
}
