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

#include "LagrangianR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "LagrangianDS.hpp"

using namespace std;

void LagrangianR::initComponents()
{
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeDS = interaction()->getSizeOfDS();

  // The initialization of Jach[0] depends on the way the Relation was built ie if the matrix
  // was read from xml or not
  if (! _jachq)
    _jachq.reset(new SimpleMatrix(sizeY, sizeDS));
  else
  {
    if (_jachq->size(0) == 0) // if the matrix dim are null
      _jachq->resize(sizeY, sizeDS);
    else
      assert((_jachq->size(1) == sizeDS && _jachq->size(0) == sizeY) &&
             "LagrangianR::initComponents inconsistent sizes between Jach[0] matrix and the interaction.");
  }
  // Added by Son Nguyen (8/12/2010)
  if (! _jachqDot)
    _jachqDot.reset(new SimpleMatrix(sizeY, sizeDS));
  else
  {
    if (_jachqDot->size(0) == 0) // if the matrix dimension are null
      _jachqDot->resize(sizeY, sizeDS);
    else
      assert((_jachqDot->size(1) == sizeDS && _jachqDot->size(0) == sizeY) &&
             "LagrangianR::initComponents inconsistent sizes between Jach[1] matrix and the interaction.");
  }
  _workX.reset(new SimpleVector(sizeDS));
  _workXdot.reset(new SimpleVector(sizeDS));
  _workZ.reset(new SimpleVector(interaction()->getSizez()));
  _workY.reset(new SimpleVector(sizeY));
}

void LagrangianR::initialize(SP::Interaction inter)
{
  assert(inter && "Lagrangian::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Memory allocation for G[i], if required (depends on the chosen constructor).
  initComponents();
  data.resize(sizeDataNames);

  LinkData();
}

void LagrangianR::LinkData()
{


  data[q0].reset(new BlockVector()); // displacement
  data[q1].reset(new BlockVector()); // velocity
  data[q2].reset(new BlockVector()); // acceleration
  data[z].reset(new BlockVector()); // z vector
  data[p0].reset(new BlockVector());
  data[p1].reset(new BlockVector());
  data[p2].reset(new BlockVector());

  SP::LagrangianDS lds;
  DSIterator it;
  for (it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    // check dynamical system type
    assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
            Type::value(**it) == Type::LagrangianDS));

    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = boost::static_pointer_cast<LagrangianDS> (*it);

    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->q());
    data[q1]->insertPtr(lds->velocity());
    data[q2]->insertPtr(lds->acceleration());
    data[z]->insertPtr(lds->z());

    // Put NonsmoothInput _p of each DS into a block. (Pointers links, no copy!!)
    for (unsigned int k = interaction()->lowerLevelForInput() ;
         k < interaction()->upperLevelForInput() + 1;
         k++)
    {
      assert(lds->p(k));
      assert(data[p0 + k]);
      data[p0 + k]->insertPtr(lds->p(k));
    }


  }
}



void LagrangianR::LinkDataFromMemory(unsigned memoryLevel)
{


  assert(memoryLevel >= 0);

  DSIterator it;
  data[q0].reset(new BlockVector()); // displacement
  data[q1].reset(new BlockVector()); // velocity
  data[q2].reset(new BlockVector()); // acceleration
  data[z].reset(new BlockVector()); // z vector
  data[p0].reset(new BlockVector());
  data[p1].reset(new BlockVector());
  data[p2].reset(new BlockVector());

  SP::LagrangianDS lds;
  for (it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    // check dynamical system type
    assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
            Type::value(**it) == Type::LagrangianDS));
    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = boost::static_pointer_cast<LagrangianDS> (*it);

    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    data[q0]->insertPtr(lds->qMemory()->getSiconosVector(memoryLevel));
    data[q1]->insertPtr(lds->velocityMemory()->getSiconosVector(memoryLevel));


    // Do nothing for the remaining of data since there are no Memory
    // An access to the content ofdata[q2] based on a link on Memory
    //must throw an exeption


    // data[q2]->insertPtr( lds->acceleration());
    // data[z]->insertPtr( lds->z());

    // // Put NonsmoothInput _p of each DS into a block. (Pointers links, no copy!!)
    // for (unsigned int k = interaction()->lowerLevelForInput() ;
    //      k < interaction()->upperLevelForInput()+1;
    //      k++)
    // {
    //   assert(lds->p(k));
    //   assert(data[p0+k]);
    //   data[p0+k]->insertPtr( lds->p(k) );
    // }

  }
}

void LagrangianR::computeh(double)
{
  RuntimeException::selfThrow("LagrangianR::computeh: not yet implemented (or useless) for Lagrangian relation of type " + subType);
}

void LagrangianR::saveRelationToXML() const
{
  RuntimeException::selfThrow("LagrangianR::saveRelationToXML - not yet implemented.");
}

SP::SiconosMatrix LagrangianR::C() const
{
  //std::cout << " SP::SiconosMatrix LagrangianR::C()      " << std::endl;
  // jachq()->display();
  //  std::cout << "jachq().get()     " <<  jachq().get()  <<std::endl;
  //return jachq();
  return _jachq;
}


void LagrangianR::display() const
{
  Relation::display();
  std::cout << " _jachq :" << std::endl;
  if (_jachq)
    _jachq->display();
  std::cout << " _jachqDot :" << std::endl;
  if (_jachqDot)
    _jachqDot->display();
  std::cout << " _jachlambda :" << std::endl;
  if (_jachlambda)
    _jachlambda->display();
  else
    std::cout << " NULL :" << std::endl;

}
