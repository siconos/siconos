/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
using namespace RELATION;

// Xml constructor
LagrangianLinearR::LagrangianLinearR(SP::RelationXML relxml):
  LagrangianR(relxml, LinearR)
{
  SP::LagrangianLinearRXML LLRxml = (boost::static_pointer_cast<LagrangianLinearRXML>(relationxml));

  // H is the minimal required input.
  if (!LLRxml->hasH())
    RuntimeException::selfThrow("LagrangianLinearR:: xml constructor failed: can not find input for H matrix.");

  H.reset(new SimpleMatrix(LLRxml->getH()));

  if (LLRxml->hasB())
  {
    b.reset(new SimpleVector(LLRxml->getB()));
  }

  if (LLRxml->hasD())
  {
    D.reset(new SimpleMatrix(LLRxml->getD()));
  }

  if (LLRxml->hasF())
  {
    F.reset(new SimpleMatrix(LLRxml->getF()));
  }
}

// Constructor from data: H, b.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB):
  LagrangianR(LinearR)
{
  H.reset(new SimpleMatrix(newH));
  b.reset(new SimpleVector(newB));
}

// Constructor from data: H.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH):
  LagrangianR(LinearR)
{
  H.reset(new SimpleMatrix(newH));
}

// Constructor from data: H, b, and D.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB, const SiconosMatrix& newD):
  LagrangianR(LinearR)
{
  H.reset(new SimpleMatrix(newH));
  D.reset(new SimpleMatrix(newD));
  b.reset(new SimpleVector(newB));
}

// Constructor from data: H, b, D and F.
LagrangianLinearR::LagrangianLinearR(const SiconosMatrix& newH, const SimpleVector& newB, const SiconosMatrix& newD, const SiconosMatrix& newF):
  LagrangianR(LinearR)
{
  H.reset(new SimpleMatrix(newH));
  D.reset(new SimpleMatrix(newD));
  b.reset(new SimpleVector(newB));
  F.reset(new SimpleMatrix(newF));
}

LagrangianLinearR::~LagrangianLinearR()
{
}

void LagrangianLinearR::initComponents()
{
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeDS = interaction->getSizeOfDS();

  assert((H->size(1) == sizeDS && H->size(0) == sizeY) &&
         "LagrangianLinearR::initComponents inconsistent sizes between H matrix and the interaction.");

  assert((if (D)
          (D->size(0) == sizeY && D->size(1) != sizeY)) &&
         "LagrangianLinearR::initComponents inconsistent sizes between D matrix and the interaction.");
  assert((if (b)
          (b->size() == sizeY)) &&
         "LagrangianLinearR::initComponents inconsistent sizes between b vector and the dimension of the interaction.");

  assert((if (F)
{
  unsigned int sizeZ = interaction->getSizeZ();
    (F->size(0) == sizeY && F->size(1) == sizeZ);
  }) &&
         "LagrangianLinearR::initComponents inconsistent sizes between F matrix and the interaction.");


  workL.reset(new SimpleVector(sizeY));

}

// Setters

void LagrangianLinearR::setH(const SiconosMatrix& newValue)
{
  if (! H)
  {
    H.reset(new SimpleMatrix(newValue));
  }
  else
    *H = newValue;
}

void LagrangianLinearR::setHPtr(SP::SiconosMatrix newPtr)
{
  H = newPtr;
}

void LagrangianLinearR::setB(const SimpleVector& newValue)
{
  if (! b)
  {
    b.reset(new SimpleVector(newValue));
  }
  else
    *b = newValue;
}

void LagrangianLinearR::setBPtr(SP::SimpleVector newPtr)
{
  b = newPtr;
}

void LagrangianLinearR::setD(const SiconosMatrix& newValue)
{
  if (! D)
    D.reset(new SimpleMatrix(newValue));
  else
    *D = newValue;
}

void LagrangianLinearR::setDPtr(SP::SiconosMatrix newPtr)
{
  D = newPtr;
}

void LagrangianLinearR::setF(const SiconosMatrix& newValue)
{
  if (! F)
  {
    F.reset(new SimpleMatrix(newValue));
  }
  else
    *F = newValue;
}

void LagrangianLinearR::setFPtr(SP::SiconosMatrix newPtr)
{
  F = newPtr;
}

void LagrangianLinearR::computeOutput(double time, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SP::SiconosVector y = interaction->getYPtr(derivativeNumber);
  SP::SiconosVector lambda = interaction->getLambdaPtr(derivativeNumber);

  //string name = "q"+toString<unsigned int>(derivativeNumber);

  if (derivativeNumber == 0)
  {
    prod(*H, *data["q0"], *y);
    if (b)
      *y += *b;
    if (D)
      prod(*D, *lambda, *y, false) ;

    if (F)
      prod(*F, *data["z"], *y, false);
  }

  else if (derivativeNumber == 1)
  {
    prod(*H, *data["q1"], *y);
    if (D)
      prod(*D, *lambda, *y, false) ;

    if (F)
      prod(*F, *data["z"], *y, false);
  }
  else if (derivativeNumber == 2)
  {
    prod(*H, *data["q2"], *y);
    if (D)
      prod(*D, *lambda, *y, false) ;

    if (F)
      prod(*F, *data["z"], *y, false);
  }

}

void LagrangianLinearR::computeInput(double time, const unsigned int level)
{
  // get lambda of the concerned interaction
  //  string name = "p"+toString<unsigned int>(level);

  *workL = *interaction->getLambdaPtr(level);
  // compute p = Ht lambda
  if (level == 0)
  {
    prod(*workL, *H, *data["p0"], false);
  }
  else if (level == 1)
  {
    prod(*workL, *H, *data["p1"], false);
  }
  else if (level == 2)
  {
    prod(*workL, *H, *data["p2"], false);
  }
  //gemv(CblasTrans,1.0,*H,*lambda,1.0, *data[name]); => not yet implemented for BlockVectors.
}

void LagrangianLinearR::saveRelationToXML() const
{
  assert(relationxml &&
         "LagrangianLinearR::saveRelationToXML - object RelationXML does not exist");

  (boost::static_pointer_cast<LagrangianLinearRXML>(relationxml))->setH(*H) ;
  (boost::static_pointer_cast<LagrangianLinearRXML>(relationxml))->setB(*b) ;
  (boost::static_pointer_cast<LagrangianLinearRXML>(relationxml))->setD(*D) ;
  (boost::static_pointer_cast<LagrangianLinearRXML>(relationxml))->setF(*F) ;
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
  if (H)
    H->display();
  else
    cout << " -> NULL " << endl;
  cout << " b: " << endl;
  if (b)
    b->display();
  else
    cout << " -> NULL " << endl;
  cout << " D: " << endl;
  if (D)
    D->display();
  else
    cout << " -> NULL " << endl;
  cout << " F: " << endl;
  if (F)
    F->display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
}
