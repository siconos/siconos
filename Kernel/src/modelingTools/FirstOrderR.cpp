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
#include "FirstOrderR.h"
#include "RelationXML.h"
#include "Interaction.h"
#include "FirstOrderNonLinearDS.h"

using namespace std;

template <class T> void FirstOrderR<T>::initDSLinks()
{
  data.resize(sizeDataNames);
  // Get the DS concerned by the interaction of this relation
  data[x].reset(new BlockVector()); // displacements
  data[z].reset(new BlockVector());
  data[r].reset(new BlockVector());

  SP::FirstOrderNonLinearDS ds;
  for (DSIterator it = interaction->dynamicalSystemsBegin(); it != interaction->dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    ds = boost::static_pointer_cast<FirstOrderNonLinearDS> (*it);
    data[x]->insertPtr(ds->getXPtr());
    data[z]->insertPtr(ds->getZPtr());
    data[r]->insertPtr(ds->getRPtr());
  }
}

template <class T> void FirstOrderR<T>::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderR::initialize failed. No Interaction linked to the present relation.");
  interaction = inter;

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = interaction->getSizeOfY();
  unsigned int sizeX = interaction->getSizeOfDS();
  unsigned int sizeZ = interaction->getSizeZ();

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  workX.reset(new SimpleVector(sizeX));
  workZ.reset(new SimpleVector(sizeZ));
  workY.reset(new SimpleVector(sizeY));
}

template <class T> void FirstOrderR<T>::setComputeHFunction(const string& pluginPath, const string& functionName)
{
  hPlugged = Plugin::setFunction(&output, pluginPath, functionName, hName);
}

template <class T> void FirstOrderR<T>::setComputeJacobianHFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  JacH[index]->setComputeFunction(pluginPath, functionName);
}

template <class T> void FirstOrderR<T>::setComputeGFunction(const string& pluginPath, const string& functionName)
{
  gPlugged = Plugin::setFunction(&input, pluginPath, functionName, gName);
}

template <class T> void FirstOrderR<T>::setComputeJacobianGFunction(const string& pluginPath, const string& functionName, unsigned int index)
{
  JacG[index]->setComputeFunction(pluginPath, functionName);
}

template <class T> void FirstOrderR<T>::computeJacH(double, unsigned int)
{
  RuntimeException::selfThrow("FirstOrderR::computeJacobianH, not (yet) implemented or forbidden for relations of type " + subType);
}

template <class T> void FirstOrderR<T>::computeJacG(double, unsigned int)
{
  RuntimeException::selfThrow("FirstOrderR::computeJacobianG, not (yet) implemented or forbidden for relations of type " + subType);
}

template <class T> void FirstOrderR<T>::display() const
{
  Relation::display();
}

