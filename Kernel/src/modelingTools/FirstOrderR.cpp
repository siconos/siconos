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
#include "FirstOrderR.hpp"
#include "RelationXML.hpp"
#include "Interaction.hpp"
#include "FirstOrderNonLinearDS.hpp"

using namespace std;

void FirstOrderR::initDSLinks()
{
  data.resize(sizeDataNames);
  // Get the DS concerned by the interaction of this relation
  data[x].reset(new BlockVector()); // displacements
  data[z].reset(new BlockVector());
  data[r].reset(new BlockVector());
  data[ds_xp].reset(new BlockVector());
  data[g_alpha].reset(new BlockVector());

  SP::FirstOrderNonLinearDS ds;
  for (DSIterator it = interaction()->dynamicalSystemsBegin(); it != interaction()->dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    ds = boost::static_pointer_cast<FirstOrderNonLinearDS> (*it);
    data[x]->insertPtr(ds->x());
    data[z]->insertPtr(ds->z());
    data[r]->insertPtr(ds->r());
    data[g_alpha]->insertPtr(ds->gAlpha());
    data[ds_xp]->insertPtr(ds->xp());
    //      data[Blambda]->insertPtr( ds->bLambda());

  }
}

void FirstOrderR::initialize(SP::Interaction inter)
{
  assert(inter && "FirstOrderR::initialize failed. No Interaction linked to the present relation.");
  _interaction = inter;

  // Check if an Interaction is connected to the Relation.
  unsigned int sizeY = interaction()->getSizeOfY();
  unsigned int sizeX = interaction()->getSizeOfDS();
  unsigned int sizeZ = interaction()->getSizez();

  // Update data member (links to DS variables)
  initDSLinks();
  // Initialize work vectors

  _workR.reset(new SimpleVector(sizeX));
  _workX.reset(new SimpleVector(sizeX));
  _workZ.reset(new SimpleVector(sizeZ));
  _workY.reset(new SimpleVector(sizeY));
}

void FirstOrderR::computeJachx(double)
{
  RuntimeException::selfThrow("FirstOrderR::computeJachx, not (yet) implemented or forbidden for relations of type " + subType);
}
void FirstOrderR::computeJachlambda(double)
{
  RuntimeException::selfThrow("FirstOrderR::computeJachlambda, not (yet) implemented or forbidden for relations of type " + subType);
}

void FirstOrderR::computeJacglambda(double t)
{
  RuntimeException::selfThrow("FirstOrderR::computeJacglambda, not (yet) implemented or forbidden for relations of type " + subType);
}

void FirstOrderR::display() const
{
  Relation::display();
}

