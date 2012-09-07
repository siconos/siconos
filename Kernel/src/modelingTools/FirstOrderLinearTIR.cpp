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
  SP::LinearRXML folrXML = cpp11ns::static_pointer_cast<LinearRXML>(_relationxml);
  // get matrices values. All are optional.

  if (folrXML->hasC())
    _jachx.reset(new SimpleMatrix(folrXML->getC()));
  else
    RuntimeException::selfThrow("FirstOrderLinearTIR:: xml constructor failed, can not find a definition for C.");

  if (folrXML->hasD())
    _jachlambda.reset(new SimpleMatrix(folrXML->getD()));

  if (folrXML->hasF())
    _F.reset(new SimpleMatrix(folrXML->getF()));

  if (folrXML->hasE())
    _e.reset(new SiconosVector(folrXML->getE()));

  if (folrXML->hasB())
    _jacglambda.reset(new SimpleMatrix(folrXML->getB()));
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
  _jachx = newC;
  _jacglambda = newB;
}

// Constructor from a complete set of data
FirstOrderLinearTIR::FirstOrderLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE, SP::SiconosMatrix newB):
  FirstOrderR(LinearTIR)
{
  _jachx = newC;
  _jacglambda = newB;
  _jachlambda = newD;
  _F = newF;
  _e = newE;
}

// Minimum data (C, B as matrices) constructor
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB):
  FirstOrderR(LinearTIR)
{
  _jachx = createSPtrSiconosMatrix((SiconosMatrix&) newC);
  _jacglambda = createSPtrSiconosMatrix((SiconosMatrix&) newB);
}


// Constructor from a complete set of data (matrices)
FirstOrderLinearTIR::FirstOrderLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE, const SiconosMatrix& newB):
  FirstOrderR(LinearTIR)
{
  _jachx = createSPtrSiconosMatrix((SiconosMatrix&) newC);
  _jacglambda = createSPtrSiconosMatrix((SiconosMatrix&) newB);
  _jachlambda = createSPtrSiconosMatrix((SiconosMatrix&) newD);
  _F = createSPtrSiconosMatrix((SiconosMatrix&) newF);
  _e.reset(new SiconosVector(newE));
}

void FirstOrderLinearTIR::initialize(Interaction & inter)
{
  // Note: do not call FirstOrderR::initialize to avoid jacobianH and jacobianG allocation.

  if (!_jachx)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() C is null and is a required input.");
  if (!_jacglambda)
    RuntimeException::selfThrow("FirstOrderLinearTIR::initialize() B is null and is a required input.");

  // Check if various operators sizes are consistent.
  // Reference: interaction.

  assert((_jachx->size(0) == inter.getSizeOfY() && _jachx->size(1) == inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and Interaction.");

  assert((_jacglambda->size(1) == inter.getSizeOfY() && _jacglambda->size(0) ==  inter.getSizeOfDS()) && "FirstOrderLinearTIR::initialize , inconsistent size between B and interaction.");

  // C and B are the minimum inputs. The others may remain null.

  if (_jachlambda)
    assert((_jachlambda->size(0) == inter.getSizeOfY() || _jachlambda->size(1) == inter.getSizeOfY()) && "FirstOrderLinearTIR::initialize , inconsistent size between C and D.");


  if (_F)
    assert(((_F->size(0) != inter.getSizeOfY()) && (_F->size(1) != inter.getSizez())) && "FirstOrderLinearTIR::initialize , inconsistent size between C and F.");
  if (_e)
    assert(_e->size() == inter.getSizeOfY() && "FirstOrderLinearTIR::initialize , inconsistent size between C and e.");
}

void FirstOrderLinearTIR::computeh(const double time, Interaction & inter)
{
  computeOutput(time, inter, 0);
}

void FirstOrderLinearTIR::computeg(const double time, Interaction & inter)
{
  computeInput(time, inter, 0);
}

void FirstOrderLinearTIR::computeOutput(const double time, Interaction & inter, unsigned int derivativeNumber)
{
  // Note that the second argument remains unamed since it is not used: for first order systems, we always compute
  // y[0]

  // We get y and lambda of the interaction (pointers)
  SiconosVector& y = *inter.y(0);
  SiconosVector& lambda = *inter.lambda(0);

  // compute y
  if (_jachx)
    prod(*_jachx, *inter.data(x), y);
  else
    y.zero();

  if (_jachlambda)
    prod(*_jachlambda, lambda, y, false);

  if (_e)
    y += *_e;

  if (_F)
    prod(*_F, *inter.data(z), y, false);
}

void FirstOrderLinearTIR::computeInput(const double time, Interaction & inter, unsigned int level)
{
  // We get lambda of the interaction (pointers)
  SiconosVector& lambda = *inter.lambda(level);
  prod(*_jacglambda, lambda, *inter.data(r), false);
}

void FirstOrderLinearTIR::display() const
{
  cout << " ===== Linear Time Invariant relation display ===== " << endl;
  cout << "| C " << endl;
  if (_jachx) _jachx->display();
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
  if (_jacglambda) _jacglambda->display();
  else cout << "->NULL" << endl;
  cout << " ================================================== " << endl;
}

void FirstOrderLinearTIR::saveRelationToXML() const
{
  //   if(!relationxml)
  RuntimeException::selfThrow("FirstOrderLinearTIR::saveRelationToXML, no yet implemented.");

  //   SP::FirstOrderLinearTIRXML folrXML = (cpp11ns::static_pointer_cast<FirstOrderLinearTIRXML>(relationxml));
  //   folrXML->setC( *_jachx );
  //   folrXML->setD( *_jachlambda );
  //   folrXML->setF( *F );
  //   folrXML->setE( *e );
  //   folrXML->setB( *_jacglambda );
}

FirstOrderLinearTIR* FirstOrderLinearTIR::convert(Relation *r)
{
  return dynamic_cast<FirstOrderLinearTIR*>(r);
}

