/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef INTERACTIONLINK_H
#define INTERACTIONLINK_H

#include "Interaction.h"
#include "SiconosConst.h"
#include "DynamicalSystem.h"
#include <iostream>
#include <vector>
#include <map>

class Interaction;
class DynamicalSystem;

/** \class InteractionLink
 *  \brief class that defines a link between 2 interactions and provides information about common DS positions
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) July 20, 2005
 *
 */

class InteractionLink
{

private:

  /* origin interaction */
  Interaction *originInteraction;

  /* connected interaction */
  Interaction *linkedInteraction;

  /* list of the common DynamicalSystem */
  std::vector<DynamicalSystem*> commonDS;

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn InteractionLink()
   *  \brief default (private) constructor
   */
  InteractionLink();

public:

  /** \fn InteractionLink((Interaction*, Interaction*, const std::vector<DynamicalSystem*>)
   *  \param Interaction *: a pointer to the origin interaction
   *  \param Interaction *: a pointer to the linked interaction
   *  \param vector<DynamicalSystem*> : list of common DS
   */
  InteractionLink(Interaction*, Interaction*, const std::vector<DynamicalSystem*>&);

  /** \fn ~InteractionLink()
   *  \brief destructor
   */
  virtual ~InteractionLink();

  // --- GETTERS/SETTERS ---

  // --- origin interaction ---

  /** \fn inline Interaction* getOriginInteractionPtr()
   *  \brief get the origin interaction pointer
   *  \return Interaction*
   */
  inline Interaction* getOriginInteractionPtr() const
  {
    return originInteraction;
  }

  /** \fn inline void setOriginInteractionPtr(Interaction* origin)
   *  \brief set the origin interaction
   *  \param Interaction*
   */
  inline void setOriginInteractionPtr(Interaction* newInter)
  {
    originInteraction = newInter;
  }

  // --- linked interaction ---

  /** \fn inline Interaction* getLinkedInteractionPtr()
   *  \brief get the linked interaction pointer
   *  \return Interaction*
   */
  inline Interaction* getLinkedInteractionPtr() const
  {
    return linkedInteraction;
  }

  /** \fn inline void setLinkedInteractionPtr(Interaction* linked)
   *  \brief set the linked interaction
   *  \param Interaction*
   */
  inline void setLinkedInteractionPtr(Interaction* newInter)
  {
    linkedInteraction = newInter;
  }

  // --- commonDS ---

  /** \fn vector<DynamicalSystem*> getCommonDS(void)
   *  \brief get the commonDS of this interactionLink
   *  \return a vector<DynamicalSystem*>
   */
  inline const std::vector<DynamicalSystem*> getCommonDS() const
  {
    return commonDS;
  }

  /** \fn void setCommonDS(const vector<DynamicalSystem*>&)
   *  \brief set the commonDS of this interactionLink
   *  \param  a vector<DynamicalSystem*>
   */
  inline void setCommonDS(const std::vector<DynamicalSystem*>& vs)
  {
    commonDS = vs;
  }

  /** \fn void display()
   * \brief interactionLink display
   */
  void display() const;

};

#endif // INTERACTIONLINK_H
