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

/*! \file InteractionManager.hpp
  \brief Definition of a class that manages dynamic interactions.
*/

#ifndef InteractionManager_h
#define InteractionManager_h

#include "Interaction.hpp"
#include "Model.hpp"
#include "DynamicalSystem.hpp"

#include "SiconosVisitor.hpp"

class InteractionManager : public SiconosVisitor
{
protected:
  void link(SP::Model model,
            SP::Interaction inter,
            SP::DynamicalSystem ds1,
            SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  void unlink(SP::Model, SP::Interaction inter);

public:
  virtual ~InteractionManager() {}

  // virtual void buildGraph(SP::Model model) = 0;
  // virtual void buildGraph(SP::DynamicalSystem body) = 0;
  // virtual void buildGraph(std::vector<SP::DynamicalSystem> bodies) = 0;
  // virtual void buildGraph(SP::SiconosContactor contactors) = 0;
  // void buildGraph(std::vector<SP::SiconosContactor> contactors);

  virtual void updateGraph() = 0;

  virtual void updateInteractions(SP::Model model) = 0;

  virtual void insertNonSmoothLaw(SP::NonSmoothLaw, int group1, int group2) = 0;
  virtual SP::NonSmoothLaw nonSmoothLaw(int group1, int group2) = 0;
};

#endif /* InteractionManager_h */
