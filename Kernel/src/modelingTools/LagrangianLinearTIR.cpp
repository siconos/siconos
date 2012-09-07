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
#include "LagrangianLinearTIR.hpp"
#include "LinearRXML.hpp"
#include "Interaction.hpp"
//
#include "LagrangianDS.hpp"

using namespace std;
using namespace RELATION;

// Xml constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::RelationXML relxml):
  LagrangianR(relxml, LinearTIR)
{
  SP::LinearRXML folrXML = cpp11ns::static_pointer_cast<LinearRXML>(relxml);
  // get matrices values. All are optional.

  if (folrXML->hasC())
    _jachq.reset(new SimpleMatrix(folrXML->getC()));
  else
    RuntimeException::selfThrow("LagrangianLinearTIR:: xml constructor failed, can not find a definition for C.");

  if (folrXML->hasD())
    _jachlambda.reset(new SimpleMatrix(folrXML->getD()));

  if (folrXML->hasF())
    _F.reset(new SimpleMatrix(folrXML->getF()));

  if (folrXML->hasE())
    _e.reset(new SiconosVector(folrXML->getE()));
}

// Minimum data (C as pointer) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC):
  LagrangianR(LinearTIR)
{
  _jachq = newC;
}

// Constructor from a complete set of data
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  _jachq = newC;
  _jachlambda = newD;
  _F = newF;
  _e = newE;
}

// Minimum data (C, e as pointers) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  _jachq = newC;
  _e = newE;
}

// Minimum data (C as matrix) constructor
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC):
  LagrangianR(LinearTIR)
{
  _jachq.reset(new SimpleMatrix(newC));
}

// Constructor from a complete set of data (matrices)
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  _jachq.reset(new SimpleMatrix(newC));
  _jachlambda.reset(new SimpleMatrix(newD));
  _F.reset(new SimpleMatrix(newF));
  _e.reset(new SiconosVector(newE));
}

// Constructor from C and e as matrix/vector
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  _jachq.reset(new SimpleMatrix(newC));
  _e.reset(new SiconosVector(newE));
}

void LagrangianLinearTIR::initComponents(Interaction& inter)
{
  unsigned int sizeY = inter.getSizeOfY();

  if (!(_jachq) || _jachq->size(1) !=  inter.getSizeOfDS() ||  _jachq->size(0) != sizeY)
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between H matrix and the interaction.");

  if ((_jachlambda) && (_jachlambda->size(0) != sizeY || _jachlambda->size(1) != sizeY))
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between D matrix and the interaction.");

  if ((_e) && _e->size() != sizeY)
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between e vector and the dimension of the interaction.");

  if ((_F) && (
        _F->size(0) != inter.getSizez() || _F->size(1) != inter.getSizez()))
    RuntimeException::selfThrow("LagrangianLinearTIR::initComponents inconsistent sizes between F matrix and the interaction.");


}

void LagrangianLinearTIR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SiconosVector& y = *inter.y(derivativeNumber);

  prod(*_jachq, *inter.data(q0 + derivativeNumber), y);

  if (derivativeNumber == 0)
  {
    if (_e)
      y += *_e;
    if (_F)
      prod(*_F, *inter.data(z), y, false);
  }

  if (_jachlambda)
  {
    SiconosVector& lambda = *inter.lambda(derivativeNumber);
    prod(*_jachlambda, lambda, y, false);
  }


}

void LagrangianLinearTIR::computeInput(const double time, Interaction& inter, unsigned int level)
{
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
  // computation of p = Ht lambda
  prod(lambda, *_jachq, *inter.data(p0 + level), false);
}

void LagrangianLinearTIR::saveRelationToXML() const
{
  assert(_relationxml &&
         "LagrangianLinearTIR::saveRelationToXML - object RelationXML does not exist");

  (cpp11ns::static_pointer_cast<LinearRXML>(_relationxml))->setC(*_jachq) ;
  (cpp11ns::static_pointer_cast<LinearRXML>(_relationxml))->setE(*_e) ;
  (cpp11ns::static_pointer_cast<LinearRXML>(_relationxml))->setD(*_jachlambda) ;
  (cpp11ns::static_pointer_cast<LinearRXML>(_relationxml))->setF(*_F) ;
}

LagrangianLinearTIR* LagrangianLinearTIR::convert(Relation *r)
{
  return dynamic_cast<LagrangianLinearTIR*>(r);
}

void LagrangianLinearTIR::display() const
{
  LagrangianR::display();
  cout << "===== Lagrangian Linear Relation display ===== " << endl;
  cout << " C: " << endl;
  if (_jachq)
    _jachq->display();
  else
    cout << " -> NULL " << endl;
  cout << " e: " << endl;
  if (_e)
    _e->display();
  else
    cout << " -> NULL " << endl;
  cout << " D: " << endl;
  if (_jachlambda)
    _jachlambda->display();
  else
    cout << " -> NULL " << endl;
  cout << " F: " << endl;
  if (_F)
    _F->display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
}
