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
