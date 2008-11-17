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
#include "FirstOrderLinearTIR.h"
#include "LinearRXML.h"
#include "Interaction.h"
#include "FirstOrderR.cpp"

using namespace std;
using namespace RELATION;

// xml constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::RelationXML relxml):
  BaseClass(relxml, LinearTIR)
{
  SP::LinearRXML folrXML = boost::static_pointer_cast<LinearRXML>(relationxml);
  // get matrices values. All are optional.

  if (folrXML->hasC())
    C.reset(new SimpleMatrix(folrXML->getC()));
  else
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor failed, can not find a definition for C.");

  if (folrXML->hasD())
    D.reset(new SimpleMatrix(folrXML->getD()));

  if (folrXML->hasF())
    F.reset(new SimpleMatrix(folrXML->getF()));

  if (folrXML->hasE())
    e.reset(new SimpleVector(folrXML->getE()));

  if (folrXML->hasB())
    B.reset(new SimpleMatrix(folrXML->getB()));
  else
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor failed, can not find a definition for B.");
}

// Minimum data (C, B as pointers) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  BaseClass(LinearTIR)
{
  C = newC;
  B = newB;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  BaseClass(LinearTIR)
{
  C = newC;
  B = newB;
  D = newD;
  F = newF;
  e = newE;
}

// Minimum data (C, B as matrices) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  BaseClass(LinearTIR)
{
  C.reset(new SimpleMatrix(newC));
  B.reset(new SimpleMatrix(newB));
}

// Constructor from a complete set of data (matrices)
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE, const SiconosMatrix& newB):
  BaseClass(LinearTIR)
{
  C.reset(new SimpleMatrix(newC));
  D.reset(new SimpleMatrix(newD));
  F.reset(new SimpleMatrix(newF));
  e.reset(new SimpleVector(newE));
  B.reset(new SimpleMatrix(newB));
}

void FirstOrderLinearTIR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderLinearTIR::initialize failed. No Interaction linked to the present relation.");
  interaction = inter;

  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Update data member (links to DS variables)
  initDSLinks();
  if (!C)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() C is null and is a required input.");
  if (!B)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeY = getInteractionPtr()->getSizeOfY();
  unsigned int sizeX = getInteractionPtr()->getSizeOfDS();
  unsigned int sizeZ = getInteractionPtr()->getSizeZ();

  assert((C->size(0) == sizeY && C->size(1) == sizeX) && "FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  assert((B->size(1) == sizeY && B->size(0) == sizeX) && "FirstOrderLinearTIR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (D)
    assert((D->size(0) == sizeY || D->size(1) == sizeY) && "FirstOrderLinearTIR::initialize , inconsistent size between C and D.");


  if (F)
    assert(((F->size(0) != sizeY) && (F->size(1) != sizeZ)) && "FirstOrderLinearTIR::initialize , inconsistent size between C and F.");
  if (e)
    assert(e->size() == sizeY && "FirstOrderLinearTIR::initialize , inconsistent size between C and e.");

  workZ.reset(new SimpleVector(sizeZ));
}

void FirstOrderLinearTIR::computeH(double time)
{
  computeOutput(time, 0);
}

void FirstOrderLinearTIR::computeG(double time)
{
  computeInput(time, 0);
}

void FirstOrderLinearTIR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = getInteractionPtr()->getYPtr(0);
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(0);

  // compute y
  if (C)
    prod(*C, *data[x], *y);
  else
    y->zero();

  if (D)
    prod(*D, *lambda, *y, false);

  if (e)
    *y += *e;

  if (F)
    prod(*F, *data[z], *y, false);
}

void FirstOrderLinearTIR::computeInput(double time, unsigned int level)
{
  // We get lambda of the interaction (pointers)
  SP::SiconosVector lambda = getInteractionPtr()->getLambdaPtr(level);
  prod(*B, *lambda, *data[r], false);
}

const SimpleMatrix FirstOrderLinearTIR::getJacH(unsigned int  index) const
{
  assert(index < 2 && "FirstOrderLinearTIR::getJacH(index) error, index is out of range.");
  if (index == 0)
    return *C;
  else
    return *D;
}

SP::SiconosMatrix FirstOrderLinearTIR::getJacHPtr(unsigned int index) const
{
  assert(index < 2 && "FirstOrderLinearTIR::getJacHPtr(index) error, index is out of range.");
  if (index == 0)
    return C;
  else
    return D;
}


const SimpleMatrix FirstOrderLinearTIR::getJacG(unsigned int  index) const
{
  assert(index < 1 && "FirstOrderLinearTIR::getJacG(index) error, index is out of range.");
  return *B;
}

SP::SiconosMatrix FirstOrderLinearTIR::getJacGPtr(unsigned int index) const
{
  assert(index < 1 && "FirstOrderLinearTIR::getJacGPtr(index) error, index is out of range.");
  return B;
}


void FirstOrderLinearTIR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (C) C->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (D) D->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (F) F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (e) e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (B) B->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearTIR::saveRelationToXML() const
{
  //   if(!relationxml)
  RuntimeException::selfThrow("FirstOrderLinearTIR::saveRelationToXML, no yet implemented.");

  //   SP::FirstOrderLinearTIRXML folrXML = (boost::static_pointer_cast<FirstOrderLinearTIRXML>(relationxml));
  //   folrXML->setC( *C );
  //   folrXML->setD( *D );
  //   folrXML->setF( *F );
  //   folrXML->setE( *e );
  //   folrXML->setB( *B );
}

FirstOrderLinearTIR* FirstOrderLinearTIR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearTIR*>(r);
}

