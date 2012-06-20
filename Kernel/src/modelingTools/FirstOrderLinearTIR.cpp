/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "FirstOrderLinearTIR.hpp"
#include "LinearRXML.hpp"
#include "Interaction.hpp"

using namespace std;
using namespace RELATION;

// xml constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::RelationXML relxml):
  FirstOrderR(relxml, LinearTIR)
{
  SP::LinearRXML folrXML = boost::static_pointer_cast<LinearRXML>(relationxml);
  // get matrices values. All are optional.

  if (folrXML->hasC())
    Jachx.reset(new SimpleMatrix(folrXML->getC()));
  else
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor failed, can not find a definition for C.");

  if (folrXML->hasD())
    _jachlambda.reset(new SimpleMatrix(folrXML->getD()));

  if (folrXML->hasF())
    _F.reset(new SimpleMatrix(folrXML->getF()));

  if (folrXML->hasE())
    _e.reset(new SiconosVector(folrXML->getE()));

  if (folrXML->hasB())
    Jacglambda.reset(new SimpleMatrix(folrXML->getB()));
  else
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor failed, can not find a definition for B.");
}
FirstOrderLinearTIR::FirstOrderLinearTIR():
  FirstOrderR(LinearTIR)
{
}
// Minimum data (C, B as pointers) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newB):
  FirstOrderR(LinearTIR)
{
  Jachx = newC;
  Jacglambda = newB;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  FirstOrderR(LinearTIR)
{
  Jachx = newC;
  Jacglambda = newB;
  _jachlambda = newD;
  _F = newF;
  _e = newE;
}

// Minimum data (C, B as matrices) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderR(LinearTIR)
{
  Jachx = createSPtrSiconosMatrix((SiconosMatrix&) newC);
  Jacglambda = createSPtrSiconosMatrix((SiconosMatrix&) newB);
}

// Constructor from a complete set of data (matrices)
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE, const SiconosMatrix& newB):
  FirstOrderR(LinearTIR)
{

  Jachx = createSPtrSiconosMatrix((SiconosMatrix&) newC);
  Jacglambda = createSPtrSiconosMatrix((SiconosMatrix&) newB);
  _jachlambda = createSPtrSiconosMatrix((SiconosMatrix&) newD);
  _F = createSPtrSiconosMatrix((SiconosMatrix&) newF);
  _e = createSPtrSiconosVector((SiconosVector&) newE);
}

void FirstOrderLinearTIR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderLinearTIR::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  // Update data member (links to DS variables)
  initDSLinks();
  if (!Jachx)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() C is null and is a required input.");
  if (!Jacglambda)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.
  unsigned int sizeZ = interaction()->getSizez();

  assert((Jachx->size(0) == interaction()->getSizeOfY() && Jachx->size(1) == interaction()->getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  assert((Jacglambda->size(1) == interaction()->getSizeOfY() && Jacglambda->size(0) ==  interaction()->getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (_jachlambda)
    assert((_jachlambda->size(0) == interaction()->getSizeOfY() || _jachlambda->size(1) == interaction()->getSizeOfY()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and D.");


  if (_F)
    assert(((_F->size(0) != interaction()->getSizeOfY()) && (_F->size(1) != sizeZ)) && "FirstOrderLinearTIR::initialize , inconsistent size between C and F.");
  if (_e)
    assert(_e->size() == interaction()->getSizeOfY() && "FirstOrderLinearTIR::initialize , inconsistent size between C and e.");

  _workZ.reset(new SiconosVector(sizeZ));
}

void FirstOrderLinearTIR::computeh(double time)
{
  computeOutput(time, 0);
}

void FirstOrderLinearTIR::computeg(double time)
{
  computeInput(time, 0);
}

void FirstOrderLinearTIR::computeOutput(double time, unsigned int)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SP::SiconosVector y = interaction()->y(0);
  SP::SiconosVector lambda = interaction()->lambda(0);

  // compute y
  if (Jachx)
    prod(*Jachx, *data[x], *y);
  else
    y->zero();

  if (_jachlambda)
    prod(*_jachlambda, *lambda, *y, false);

  if (_e)
    *y += *_e;

  if (_F)
    prod(*_F, *data[z], *y, false);
}

void FirstOrderLinearTIR::computeInput(double time, unsigned int level)
{
  // We get lambda of the interaction (pointers)
  SP::SiconosVector lambda = interaction()->lambda(level);
  prod(*Jacglambda, *lambda, *data[r], false);
}

void FirstOrderLinearTIR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (Jachx) Jachx->display();
  else cout << "->NULL" << endl;
  cout << "| D " << endl;
  if (_jachlambda) _jachlambda->display();
  else cout << "->NULL" << endl;
  cout << "| F " << endl;
  if (_F) _F->display();
  else cout << "->NULL" << endl;
  cout << "| e " << endl;
  if (_e) _e->display();
  else cout << "->NULL" << endl;
  cout << "| B " << endl;
  if (Jacglambda) Jacglambda->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearTIR::saveRelationToXML() const
{
  //   if(!relationxml)
  RuntimeException::selfThrow("FirstOrderLinearTIR::saveRelationToXML, no yet implemented.");

  //   SP::FirstOrderLinearTIRXML folrXML = (boost::static_pointer_cast<FirstOrderLinearTIRXML>(relationxml));
  //   folrXML->setC( *Jachx );
  //   folrXML->setD( *_jachlambda );
  //   folrXML->setF( *F );
  //   folrXML->setE( *e );
  //   folrXML->setB( *Jacglambda );
}

FirstOrderLinearTIR* FirstOrderLinearTIR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearTIR*>(r);
}

