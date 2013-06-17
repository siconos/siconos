/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "NonSmoothDynamicalSystem.hpp"
#include "NonSmoothDynamicalSystemXML.hpp"
#include "DynamicalSystemXML.hpp"
#include "Topology.hpp"
#include "Interaction.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "InteractionXML.hpp"


using namespace RELATION;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(): _BVP(false), _mIsLinear(true)
{
  // === Builds an empty topology ===
  _topology.reset(new Topology());
};


// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(SP::NonSmoothDynamicalSystemXML newNsdsxml):
  _BVP(false), _nsdsxml(newNsdsxml), _mIsLinear(true)
{
  assert(_nsdsxml && "NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");

  // === DS Vector fill-in: we sweep the list of DSXML and for each of
  // them, add a new DynamicalSystem in the set allDS. ===
  SetOfDSXML dsList = _nsdsxml->getDynamicalSystemsXML();
  SetOfDSXML::iterator it;
  Type::Siconos type;
  // === Builds an empty topology ===
  _topology.reset(new Topology());

  // Add dynamical systems into the nsds (i.e. into the topology and its graphs).
  for (it = dsList.begin(); it != dsList.end(); ++it)
  {
    type = (*it)->getType();
    if (type  == Type::LagrangianDS)  // LagrangianDS
      insertDynamicalSystem(SP::LagrangianDS(new LagrangianDS(*it)));
    else if (type == Type::LagrangianLinearTIDS)  // Lagrangian Linear Time Invariant
      insertDynamicalSystem(SP::LagrangianLinearTIDS(new LagrangianLinearTIDS(*it)));
    else if (type == Type::FirstOrderLinearDS)  // Linear DS
      insertDynamicalSystem(SP::FirstOrderLinearDS(new FirstOrderLinearDS(*it)));
    else if (type == Type::FirstOrderLinearTIDS)  // Linear Time Invariant DS
      insertDynamicalSystem(SP::FirstOrderLinearTIDS(new FirstOrderLinearTIDS(*it)));
    else if (type == Type::FirstOrderNonLinearDS)  // Non linear DS
      insertDynamicalSystem(SP::FirstOrderNonLinearDS(new FirstOrderNonLinearDS(*it)));
    else RuntimeException::
      selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + type);
  }
  
  // ===  The same process is applied for Interactions ===
  SetOfInteractionsXML  interactionsList = _nsdsxml->getInteractionsXML();
  SetOfInteractionsXMLIt it2;
  CheckInsertInteraction checkInter;

  for (it2 = interactionsList.begin(); it2 != interactionsList.end(); ++it2)
  {
    // Create an interaction
    SP::Interaction inter(new Interaction(*it2));
    SP::DynamicalSystem ds1, ds2;
    std::vector<int> dsNumbers;
    // Get the numbers of ds involved in the interaction (max == 2)
    (*it2)->getDSNumbers(dsNumbers);
    ds1 = dynamicalSystem(dsNumbers[0]);
    if(dsNumbers.size() == 2)
    {
      ds2 =  dynamicalSystem(dsNumbers[1]);
      link(inter, ds1, ds2);
    }
    else
      link(inter,ds1);
  }
}

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  clear();
}

// === DynamicalSystems management ===

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  if (_nsdsxml)
  {
    _nsdsxml->setBVP(_BVP);// no need to change the value of BVP, it mustn't change anyway

    DynamicalSystemsGraph::VIterator vi;
    for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
    {
      SP::DynamicalSystem ds = dynamicalSystems()->bundle(*vi);
      if (Type::value(*ds) == Type::LagrangianDS)
        (std11::static_pointer_cast<LagrangianDS>(ds))->saveDSToXML();
      else if (Type::value(*ds) == Type::LagrangianLinearTIDS)
        (std11::static_pointer_cast<LagrangianLinearTIDS>(ds))->saveDSToXML();
      else if (Type::value(*ds) == Type::FirstOrderLinearDS)
        (std11::static_pointer_cast<FirstOrderLinearDS>(ds))->saveDSToXML();
      else if (Type::value(*ds) == Type::FirstOrderLinearTIDS)
        (std11::static_pointer_cast<FirstOrderLinearDS>(ds))->saveDSToXML();
      else if (Type::value(*ds) == Type::FirstOrderNonLinearDS)
        ds->saveDSToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of DynamicalSystem");
    }

    InteractionsIterator it2;
    for (it2 = _topology->interactions()->begin();
         it2 != interactions()->end(); ++it2)
      (*it2)->saveInteractionToXML();
  }
  else RuntimeException::
    selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NonSmoothDynamicalSystemXML object doesn't exists");
}

void NonSmoothDynamicalSystem::display() const
{
  std::cout << " ===== Non Smooth Dynamical System display ===== " <<std::endl;
  std::cout << "---> isBVP = " << _BVP <<std::endl;
  dynamicalSystems()->begin();
  _topology->interactions()->display();
  std::cout << "===================================================" <<std::endl;
}

#include <limits>
double NonSmoothDynamicalSystem::nsdsConvergenceIndicator()
{
  // calculate the max value of all DS convergence indicators
  double convergenceIndicator = -std::numeric_limits<double>::infinity();
  double dsIndic ;
  DynamicalSystemsGraph::VIterator vi;
  for (vi = dynamicalSystems()->begin(); vi != dynamicalSystems()->end(); ++vi)
  {
    dsIndic = dynamicalSystems()->bundle(*vi)->dsConvergenceIndicator();
    if (dsIndic > convergenceIndicator) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}

void NonSmoothDynamicalSystem::clear()
{
  _topology->clear();
}

void NonSmoothDynamicalSystem::setControlInput(SP::DynamicalSystem ds, SP::SiconosMatrix B)
{
  DynamicalSystemsGraph& DSG0 = *_topology->dSG(0);
  DSG0.B[DSG0.descriptor(ds)] = B;
}

void NonSmoothDynamicalSystem::setObserverInput(SP::DynamicalSystem ds, SP::SiconosMatrix L)
{
  DynamicalSystemsGraph& DSG0 = *_topology->dSG(0);
  DSG0.L[DSG0.descriptor(ds)] = L;
}
