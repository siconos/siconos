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
    F.reset(new SimpleMatrix(folrXML->getF()));

  if(folrXML->hasE())
  e.reset(new SimpleVector(folrXML->getE()));*/
}

// Minimum data (C as pointer) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC):
  LagrangianR(LinearTIR)
{
  Jacqh = newC;
}

// Constructor from a complete set of data
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosMatrix newD, SP::SiconosMatrix newF, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  Jacqh = newC;
  JacLH = newD;
  F = newF;
  e = newE;
}

// Minimum data (C, e as pointers) constructor
LagrangianLinearTIR::LagrangianLinearTIR(SP::SiconosMatrix newC, SP::SiconosVector newE):
  LagrangianR(LinearTIR)
{
  Jacqh = newC;
  e = newE;
}

// Minimum data (C as matrix) constructor
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC):
  LagrangianR(LinearTIR)
{
  Jacqh.reset(new SimpleMatrix(newC));
}

// Constructor from a complete set of data (matrices)
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD, const SiconosMatrix& newF, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  Jacqh.reset(new SimpleMatrix(newC));
  JacLH.reset(new SimpleMatrix(newD));
  F.reset(new SimpleMatrix(newF));
  e.reset(new SimpleVector(newE));
}

// Constructor from C and e as matrix/vector
LagrangianLinearTIR::LagrangianLinearTIR(const SiconosMatrix& newC, const SiconosVector& newE):
  LagrangianR(LinearTIR)
{
  RuntimeException::selfThrow("LagrangianLinearTIR::LagrangianLinearTIR,  copy matrix in constructor\n");
  Jacqh.reset(new SimpleMatrix(newC));
  e.reset(new SimpleVector(newE));
}

void LagrangianLinearTIR::initComponents()
{
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();

  assert((Jacqh) ? (Jacqh->size(1) == sizeDS && Jacqh->size(0) == sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between H matrix and the interaction.");

  assert((JacLH) ? (JacLH->size(0) == sizeY && JacLH->size(1) != sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between D matrix and the interaction.");
  assert((e) ? (e->size() == sizeY) : 1 &&
         "LagrangianLinearTIR::initComponents inconsistent sizes between e vector and the dimension of the interaction.");

  assert((F) ?
         (F->size(0) == interaction()->getSizeZ() && F->size(1) == interaction()->getSizeZ()) : 1 &&
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
  prod(*Jacqh, *data[q0 + derivativeNumber], *y);

  if (derivativeNumber == 0)
  {
    if (e)
      *y += *e;
    if (F)
      prod(*F, *data[z], *y, false);
  }

  if (JacLH)
    prod(*JacLH, *lambda, *y, false) ;


}

void LagrangianLinearTIR::computeInput(double time, const unsigned int level)
{
  // get lambda of the concerned interaction
  //  string name = "p"+toString<unsigned int>(level);

  *_workL = *interaction()->lambda(level);
  // computation of p = Ht lambda
  prod(*_workL, *Jacqh, *data[p0 + level], false);
  //gemv(CblasTrans,1.0,*H,*lambda,1.0, *data[name]); => not yet implemented for BlockVectors.
}
/*
const SimpleMatrix LagrangianLinearTIR::getJacXH(unsigned int  index ) const
{
  assert(index<2&&"LagrangianLinearTIR::getJacH(index) error, index is out of range.");
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

  (boost::static_pointer_cast<LinearRXML>(relationxml))->setC(*Jacqh) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setE(*e) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setD(*JacLH) ;
  (boost::static_pointer_cast<LinearRXML>(relationxml))->setF(*F) ;
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
  if (Jacqh)
    Jacqh->display();
  else
    cout << " -> NULL " << endl;
  cout << " e: " << endl;
  if (e)
    e->display();
  else
    cout << " -> NULL " << endl;
  cout << " D: " << endl;
  if (JacLH)
    JacLH->display();
  else
    cout << " -> NULL " << endl;
  cout << " F: " << endl;
  if (F)
    F->display();
  else
    cout << " -> NULL " << endl;
  cout << "===================================== " << endl;
}
