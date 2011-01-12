/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

void NewtonEulerR::initComponents()
{
  _ysize = interaction()->getSizeOfY();
  _xsize = interaction()->getSizeOfDS();
  _qsize = 7 * (_xsize / 6);

  // The initialization of Jach[0] depends on the way the Relation was built ie if the matrix
  // was read from xml or not
  if (! _jachq)
    _jachq.reset(new SimpleMatrix(_ysize, _qsize));
  else
  {
    if (_jachq->size(0) == 0) // if the matrix dim are null
    {
      _jachq->resize(_ysize, _qsize);
    }
    else
    {
      assert((_jachq->size(1) == _qsize && _jachq->size(0) == _ysize) ||
             (printf("NewtonEuler::initComponents _jachq->size(1) = %d ,_qsize = %d , _jachq->size(0) = %d ,_ysize =%d \n", _jachq->size(1), _qsize, _jachq->size(0), _ysize) && false) ||
             ("NewtonEuler::initComponents inconsistent sizes between _jachq matrix and the interaction." && false));
    }
  }
  if (! _jachqT)
    _jachqT.reset(new SimpleMatrix(_ysize, _xsize));


  _workX.reset(new SimpleVector(_xsize));
  _workQ.reset(new SimpleVector(_qsize));
  _workZ.reset(new SimpleVector(interaction()->getSizez()));
  _workY.reset(new SimpleVector(_ysize));
  //_yQ.reset(new SimpleVector(1));
  _jachqProj = _jachq;
  SP::BlockVector vaux = interaction()->y(0);
  SP::SiconosVector vaux2 = (*vaux)[0];
  _yProj = boost::static_pointer_cast<SimpleVector>((*(interaction()->y(0)))[0]);
}

void NewtonEulerR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderR::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents();
  data.resize(sizeDataNames);

  DSIterator it;
  data[q0].reset(new BlockVector()); // displacement
  data[deltaq].reset(new BlockVector());
  data[q1].reset(new BlockVector()); // velocity
  //  data[q2].reset(new BlockVector()); // acceleration
  data[z].reset(new BlockVector()); // z vector
  data[p0].reset(new BlockVector());
  data[p1].reset(new BlockVector());
  data[p2].reset(new BlockVector());
  SP::NewtonEulerDS lds;
  int sizeForAllxInDs = 0;
  for (it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    Type::Siconos type = Type::value(**it);
    // check dynamical system type
    assert((type == Type::NewtonEulerDS) && "NewtonEulerR::initialize failed, not implemented for dynamical system of type: " + type);

    // convert vDS systems into NewtonEulerDS and put them in vLDS
    lds = boost::static_pointer_cast<NewtonEulerDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->q());
    data[deltaq]->insertPtr(lds->deltaq());
    data[q1]->insertPtr(lds->dotq());
    //    data[q2]->insertPtr( lds->acceleration());
    data[p0]->insertPtr(lds->p(1));
    data[p1]->insertPtr(lds->p(1));
    data[p2]->insertPtr(lds->p(2));
    data[z]->insertPtr(lds->z());
    sizeForAllxInDs += lds->p(1)->size();
  }
  _contactForce.reset(new SimpleVector(sizeForAllxInDs));
  _contactForce->zero();
}


void NewtonEulerR::computeh(double)
{
  SP::SiconosVector y = interaction()->y(0);
  *_workQ = *data[q0];
  //prod(*_jachq,*data[q0],*y);

  prod(*_jachq, *_workQ, *y, true);
  if (_e)
    *y += *_e;
  //  printf("NewtonEulerR::computeh() q:\n");
  //  _workQ->display();
  //  printf("NewtonEulerR::computeh() :\n");
  //  y->display();
}

//  void NewtonEulerR::computeJachx(double)
// {
//   RuntimeException::selfThrow("FirstOrderR::computeJacobianXH, not (yet) implemented or forbidden for relations of type "+subType);
// }
//  void NewtonEulerR::computeJachlambda(double)
// {
//   RuntimeException::selfThrow("FirstOrderR::computeJacobianLH, not (yet) implemented or forbidden for relations of type "+subType);
// }

void NewtonEulerR::saveRelationToXML() const
{
  RuntimeException::selfThrow("NewtonEulerR1::saveRelationToXML - not yet implemented.");
}

void NewtonEulerR::display() const
{
  Relation::display();
}

// void NewtonEulerR::setqProjectOnConstraints(SP::SiconosVector lambda){
//   SP::SimpleMatrix aux(new SimpleMatrix(*_jachq));
//   aux->trans();
//   prod(*aux,*lambda,*data[q0],false);
// }

void NewtonEulerR::computeOutput(double t, unsigned int derivativeNumber)
{

  if (derivativeNumber == 0)
  {
    computeh(t);
  }
  else
  {
    computeJachq(t);

    SP::SiconosVector y = interaction()->y(derivativeNumber);
    if (derivativeNumber == 1)
    {
      prod(*_jachq, *data[q1], *y);
      //_jachq->display();
      //y->display();
      /*  if (false){
        unsigned int k=0;
        DSIterator itDS;
        itDS=interaction()->dynamicalSystemsBegin();
        while(  itDS!=interaction()->dynamicalSystemsEnd() ){
          SP::NewtonEulerDS d =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
          SP::SimpleVector vPredictor = d->vPredictor();
          unsigned int d_dim=d->getDim();
          SP::SimpleVector F(new SimpleVector(d_dim));
          setBlock(_contactForce,F,d_dim,0,k);
          d->luM()->PLUForwardBackwardInPlace(*F);
          *(d->vPredictor())+=*F;
          k+=(*itDS)->getDim();
        }
        }*/
    }
    else  //if(derivativeNumber == 2)
      //  prod(*_jachq,*data[q2],*y); // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
      //   else
      RuntimeException::selfThrow("NewtonEulerR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

/** to compute p
 *  \param double : current time
 *  \param unsigned int: "derivative" order of lambda used to compute input
 */
void NewtonEulerR::computeInput(double t, unsigned int level)
{
  /*implemented for the bouncing ball*/


  // computeJachq(t);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction()->lambda(level);

  prod(*lambda, *_jachqT, *_contactForce, true);
  /*data is a pointer of memory associated to a dynamical system*/
  /** false because it consists in doing a sum*/
  prod(*lambda, *_jachqT, *data[p0 + level], false);
}
void NewtonEulerR::computeJachqT()
{
  unsigned int k = 0;
  DSIterator itDS;
  int sizey = interaction()->getSizeOfY();
  SP::SimpleMatrix auxBloc(new SimpleMatrix(sizey, 7));
  SP::SimpleMatrix auxBloc2(new SimpleMatrix(sizey, 6));
  Index dimIndex(2);
  Index startIndex(4);
  itDS = interaction()->dynamicalSystemsBegin();
  while (itDS != interaction()->dynamicalSystemsEnd())
  {
    startIndex[0] = 0;
    startIndex[1] = 7 * k / 6;
    startIndex[2] = 0;
    startIndex[3] = 0;
    dimIndex[0] = sizey;
    dimIndex[1] = 7;
    setBlock(_jachq, auxBloc, dimIndex, startIndex);
    SP::NewtonEulerDS d =  boost::static_pointer_cast<NewtonEulerDS> (*itDS);
    SP::SiconosMatrix T = d->T();

    prod(*auxBloc, *T, *auxBloc2);

    startIndex[0] = 0;
    startIndex[1] = 0;
    startIndex[2] = 0;
    startIndex[3] = k;
    dimIndex[0] = sizey;
    dimIndex[1] = 6;

    setBlock(auxBloc2, _jachqT, dimIndex, startIndex);
    k += (*itDS)->getDim();
    itDS++;
  }
}
