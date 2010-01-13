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
      assert((_jachq->size(1) == _qsize && _jachq->size(0) == _ysize) &&
             "NewtonEuler::initComponents inconsistent sizes between Jach[0] matrix and the interaction.");
  }
  if (! _jachqT)
    _jachqT.reset(new SimpleMatrix(_ysize, _xsize));


  _workX.reset(new SimpleVector(_xsize));
  _workQ.reset(new SimpleVector(_qsize));
  _workZ.reset(new SimpleVector(interaction()->getSizeZ()));
  _workY.reset(new SimpleVector(_ysize));
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
  data[q1].reset(new BlockVector()); // velocity
  //  data[q2].reset(new BlockVector()); // acceleration
  data[z].reset(new BlockVector()); // z vector
  data[p0].reset(new BlockVector());
  data[p1].reset(new BlockVector());
  data[p2].reset(new BlockVector());
  SP::NewtonEulerDS lds;
  DS::TYPES type;
  for (it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    type = (*it)->getType();
    // check dynamical system type
    assert((type == DS::NENLDS) && "NewtonEulerR::initialize failed, not implemented for dynamical system of type: " + type);

    // convert vDS systems into NewtonEulerDS and put them in vLDS
    lds = boost::static_pointer_cast<NewtonEulerDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->q());
    data[q1]->insertPtr(lds->dotq());
    //    data[q2]->insertPtr( lds->acceleration());
    data[p0]->insertPtr(lds->p(1));
    data[p1]->insertPtr(lds->p(1));
    data[p2]->insertPtr(lds->p(2));
    data[z]->insertPtr(lds->z());
  }
}


void NewtonEulerR::computeh(double)
{
  SP::SiconosVector y = interaction()->y(0);
  *_workQ = *data[q0];
  //prod(*_jachq,*data[q0],*y);
  prod(*_jachq, *_workQ, *y);
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

void NewtonEulerR::computeOutput(double t, unsigned int derivativeNumber)
{
  /*implemented for the bouncing ball*/
  if (derivativeNumber == 0)
  {
    computeh(t);
  }
  else
  {
    SP::SiconosVector y = interaction()->y(derivativeNumber);
    if (derivativeNumber == 1)
      prod(*_jachq, *data[q1], *y);
    else //if(derivativeNumber == 2)
      //  prod(*_jachq,*data[q2],*y); // Approx: y[2] = Jach[0]q[2], other terms are neglected ...
      //   else
      RuntimeException::selfThrow("LagrangianCompliantR::computeOutput(time,index), index out of range or not yet implemented.");
  }
}

/** to compute p
 *  \param double : current time
 *  \param unsigned int: "derivative" order of lambda used to compute input
 */
void NewtonEulerR::computeInput(double t, unsigned int level)
{
  /*implemented for the bouncing ball*/



  //  computeJachq(time);
  // get lambda of the concerned interaction
  SP::SiconosVector lambda = interaction()->lambda(level);

  // data[name] += trans(G) * lambda
  prod(*lambda, *_jachqT, *data[p0 + level], false);
}
