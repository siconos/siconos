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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "NewtonEulerR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "NewtonEulerDS.hpp"

using namespace std;

//#define NER_DEBUG


void NewtonEulerR::initComponents(Interaction& inter)
{
  unsigned int ySize = inter.getSizeOfY();
  unsigned int xSize = inter.getSizeOfDS();
  unsigned int qSize = 7 * (xSize / 6);

  // The initialization of Jach[0] depends on the way the Relation was built ie if the matrix
  // was read from xml or not


  if (! _jachq)
    _jachq.reset(new SimpleMatrix(ySize, qSize));
  else
  {
    if (_jachq->size(0) == 0)
    {
      // if the matrix dim are null
      _jachq->resize(ySize, qSize);
    }
    else
    {
      assert((_jachq->size(1) == qSize && _jachq->size(0) == ySize) ||
             (printf("NewtonEuler::initComponents _jachq->size(1) = %d ,_qsize = %d , _jachq->size(0) = %d ,_ysize =%d \n", _jachq->size(1), qSize, _jachq->size(0), ySize) && false) ||
             ("NewtonEuler::initComponents inconsistent sizes between _jachq matrix and the interaction." && false));
    }
  }
#ifdef NER_DEBUG
  std::cout << "NewtonEulerR::initComponents() _jachq" << std::endl;
  _jachq->display();
#endif

  if (! _jachqT)
    _jachqT.reset(new SimpleMatrix(ySize, xSize));

}

void NewtonEulerR::initialize(Interaction& inter)
{
  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents(inter);

  _contactForce.reset(new SiconosVector(inter.data(p1)->size()));
  _contactForce->zero();
}


void NewtonEulerR::computeh(const double time, Interaction& inter)
{
  SiconosVector& y = *inter.y(0);

  prod(*_jachq, *inter.data(q0), y, true);
  if (_e)
    y += *_e;
}


void NewtonEulerR::saveRelationToXML() const
{
  RuntimeException::selfThrow("NewtonEulerR1::saveRelationToXML - not yet implemented.");
}

void NewtonEulerR::computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber)
{

  if (derivativeNumber == 0)
  {
    computeh(time, inter);
  }
  else
  {
    computeJachq(time, inter);

    SiconosVector& y = *inter.y(derivativeNumber);
    if (derivativeNumber == 1)
    {
      prod(*_jachqT, *inter.data(velo), y);
    }
    else //if(derivativeNumber == 2)
      RuntimeException::selfThrow("NewtonEulerR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

/** to compute p
*  \param double : current time
*  \Param unsigned int: "derivative" order of lambda used to compute input
*/
void NewtonEulerR::computeInput(const double time, Interaction& inter, unsigned int level)
{
  /*implemented for the bouncing ball*/


  // computeJachq(t);
  // get lambda of the concerned interaction
  SiconosVector& lambda = *inter.lambda(level);
#ifdef NER_DEBUG
  printf("\n");
  printf("NewtonEulerR::computeInput start for level %i:", level);
  std::cout << "lambda( "  << level << ")" << std::endl;
  lambda->display();
  std::cout << "data[p0+level] before " << std::endl;
  inter.data(p0 + level)->display();
#endif
  if (level == 1) /* \warning : we assume that ContactForce is given by lambda[1] */
  {
    prod(lambda, *_jachqT, *_contactForce, true);
#ifdef NER_DEBUG
    printf("NewtonEulerR::computeInput contact force :");
    _contactForce->display();
#endif

    /*data is a pointer of memory associated to a dynamical system*/
    /** false because it consists in doing a sum*/
    prod(lambda, *_jachqT, *inter.data(p0 + level), false);
  }
  else if (level == 0)
  {
    prod(lambda, *_jachq, *inter.data(p0 + level), false);
#ifdef NER_DEBUG
    std::cout << "_jachq" << std::endl;
    _jachq->display();
    std::cout << "data[p0+level]" << inter.data(p0 + level) <<  std::endl;
    std::cout << "data[p0+level]->vector(0)" << inter.data(p0 + level)->vector(0) <<  std::endl;
    if (inter.data(p0 + level)->getNumberOfBlocks() > 1)
      std::cout << "data[p0+level]->vector(1)" << inter.data(p0 + level)->vector(1) <<  std::endl;
    inter.data(p0 + level)->display();


    SP::SiconosVector buffer(new SiconosVector(inter.data(p0 + level)->size()));
    prod(*lambda, *_jachq, *buffer, true);
    std::cout << "added part to p" << buffer <<  std::endl;
    buffer->display();

    printf("NewtonEulerR::computeInput end for level %i:", level);
    printf("\n");
#endif

  }
  else
    RuntimeException::selfThrow("NewtonEulerR::computeInput(double t, unsigned int level) - not yet implemented for level > 1");
}
/*It computes _jachqT=_jachq*T. Uploaded in the case of an unilateral constraint (NewtonEulerFrom3DLocalFrameR and NewtonEulerFrom1DLocalFrameR)*/
void NewtonEulerR::computeJachqT(Interaction& inter)
{
  unsigned int k = 0;
  DSIterator itDS;
  unsigned int ySize = inter.getSizeOfY();
  SP::SimpleMatrix auxBloc(new SimpleMatrix(ySize, 7));
  SP::SimpleMatrix auxBloc2(new SimpleMatrix(ySize, 6));
  Index dimIndex(2);
  Index startIndex(4);
  itDS = inter.dynamicalSystemsBegin();
  while (itDS != inter.dynamicalSystemsEnd())
  {
    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = ySize;
    dimIndex[1] = 7;
    setBlock(_jachq, auxBloc, dimIndex, startIndex);
    NewtonEulerDS& d = *cpp11ns::static_pointer_cast<NewtonEulerDS> (*itDS);
    SiconosMatrix& T = *d.T();

    prod(*auxBloc, T, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = ySize;
    dimIndex[1] = 6;

    setBlock(auxBloc2, _jachqT, dimIndex, startIndex);
    k += (*itDS)->getDim();
    itDS++;
  }
}
