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
  /*  SP::LinearRXML folrXML = boost::static_pointer_cast<LinearRXML>(relationxml);
  // get matrices values. All are optional.

  if(folrXML->hasC())
    C.reset(new SimpleMatrix(folrXML->getC()));
  else
    RuntimeException::selfThrow("LagrangianLinearTIR:: xml constructor failed, can not find a definition for C.");

  if(folrXML->hasD())
    D.reset(new SimpleMatrix(folrXML->getD()));

  if(folrXML->hasF())
    _F.reset(new SimpleMatrix(folrXML->getF()));

  if(folrXML->hasE())
  e.reset(new SimpleVector(folrXML->getE()));*/
}

// Minimum data (C as pointer) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC):
  LagrangianR(LinearTIR)
{
  Jachq = newC;
}

// Constructor from a complete set of data
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  Jachq = newC;
  Jachlambda = newD;
  _F = newF;
  _e = newE;
}

// Minimum data (C, e as pointers) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  Jachq = newC;
  _e = newE;
}

// Minimum data (C as matrix) constructor
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC):
  LagrangianR(LinearTIR)
{
  Jachq.reset(new SimpleMatrix(newC));
}

// Constructor from a complete set of data (matrices)
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  Jachq.reset(new SimpleMatrix(newC));
  Jachlambda.reset(new SimpleMatrix(newD));
  _F.reset(new SimpleMatrix(newF));
  _e.reset(new SimpleVector(newE));
}

// Constructor from C and e as matrix/vector
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  Jachq.reset(new SimpleMatrix(newC));
  _e.reset(new SimpleVector(newE));
}

void LagrangianLinearTIR::initComponents()
{
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();

  assert((Jachq) ? (Jachq->size(1) == sizeDS && Jachq->size(0) == sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between H matrix and the interaction.");

  assert((Jachlambda) ? (Jachlambda->size(0) == sizeY && Jachlambda->size(1) != sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between D matrix and the interaction.");
  assert((_e) ? (_e->size() == sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between e vector and the dimension of the interaction.");

  assert((_F) ?
         (_F->size(0) == interaction()->getSizeZ() && _F->size(1) == interaction()->getSizeZ()) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between F matrix and the interaction.");


  _workL.reset(new SimpleVector(sizeY));

}

void LagrangianLinearTIR::computeh(double time)
{
  computeOutput(time, 0);
}
void LagrangianLinearTIR::computeg(double time)
{
  computeInput(time, 0);
}

void LagrangianLinearTIR::computeOutput(double time, unsigned int derivativeNumber)
{
  // get y and lambda of the interaction
  SP::SiconosVector y = interaction()->y(derivativeNumber);
  SP::SiconosVector lambda = interaction()->lambda(derivativeNumber);

  //string name = "q"+toString<unsigned int>(derivativeNumber);
  prod(*Jachq, *data[q0 + derivativeNumber], *y);

  if (derivativeNumber == 0)
  {
    if (_e)
      *y += *_e;
    if (_F)
      prod(*_F, *data[z], *y, false);
  }

  if (Jachlambda)
    prod(*Jachlambda, *lambda, *y, false) ;


}

void LagrangianLinearTIR::computeInput(double time, const unsigned int level)
{
  // get lambda of the concerned interaction
  //  string name = "p"+toString<unsigned int>(level);

  *_workL = *interaction()->lambda(level);
  // computation of p = Ht lambda
  prod(*_workL, *Jachq, *data[p0 + level], false);
  //gemv(CblasTrans,1.0,*H,*lambda,1.0, *data[name]); => not yet implemented for BlockVectors.
}
/*
const SimpleMatrix LagrangianLinearTIR::getJachx(unsigned int  index ) const
{
  assert(index<2&&"LagrangianLinearTIR::getJach(index) error, index is out of range.");
  if(index == 0)
    return *C;
  else
    return *D;
}
*/


void LagrangianLinearTIR::saveRelationToXML() const
{
  assert(relationxml &&
         "LagrangianLinearTIR::saveRelationToXML - object RelationXML does not exist");

  (boost::static_pointer_cast<LinearRXML>(relationxml))->setC(*Jachq) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setE(*_e) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setD(*Jachlambda) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setF(*_F) ;
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
  if (Jachq)
    Jachq->display();
  else
    cout << " -> NULL " << endl;
  cout << " e: " << endl;
  if (_e)
    _e->display();
  else
    cout << " -> NULL " << endl;
  cout << " D: " << endl;
  if (Jachlambda)
    Jachlambda->display();
  else
    cout << " -> NULL " << endl;
  cout << " F: " << endl;
  if (_F)
    _F->display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
}
