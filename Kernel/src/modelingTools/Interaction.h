/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/*! \file Interaction.h
  Interaction class and related typedef
*/


#ifndef INTERACTION_H
#define INTERACTION_H

// for composition ...
#include "NonSmoothLaw.h"
#include "Relation.h"
#include "NonSmoothDynamicalSystem.h"
#include "DynamicalSystemsSet.h"

// IO (XML)
#include "InteractionXML.h"

// const
#include "SiconosConst.h"

// stl tools
#include <vector>

class DynamicalSystemsSet;
class NonSmoothLaw;
class DynamicalSystem;
class Relation;
class NonSmoothDynamicalSystem;
class InteractionXML;


/** container for SiconosVectors */
typedef std::vector< BlockVector* > VectorOfBlocks;

/** iterator through vector of SimpleVector */
typedef VectorOfBlocks::iterator VectorOfBlocksIterator;

/** type used for inside-class allocation checking */
typedef std::deque<bool>  AllocationFlags;

//!  An Interaction describes the non-smooth interactions between a set of Dynamical Systems.
/**
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) Apr 29, 2004
 *
 * An interaction represents the "link" between a set of Dynamical Systems (var: involvedDS) that interact alltogether through
 * some relations (between state variables (x,R) and local variables (y,lambda)) completed by a non-smooth law.
 *
 * Thus, the interaction main members are:
 *
 * - a set of Dynamical Systems (from 1 to ...) that interacts, named involvedDS.
 *
 * - relation: a pointer to a Relation object that determines the type of relation and so the way it is computed.
 *   Warning: there is only one Relation object (ie only one type of relation for an interaction) but there can be several "relations", in the sense
 *   of constraints equations between (y,lambda) and (x,r).
 *
 * - nslaw: the non smooth law
 *
 * - the local variables y and lambda (their size is interactionSize).
 *   stl vectors are used and y[i] (resp lambda[i]) represents the i-eme derivative of variable y (resp lambda).
 *
 *   y is a container of BlockVector. Each block corresponds to a "unitary relation", a relation which size is the one of the non-smooth law.
 *    => ySize = interactionSize = numberOfRelations * nsLawSize .
 *   Same thing for lambda.
 *
 *  => all the relations of the interaction have the same non-smooth law.
 *  => the number of relations is equal to interactionSize/nsLawSize. Thus a relation is not necessarily represented
 *     by a single equation.
 *
 * \todo save a number of derivatives of y depending on the relative degree (at the time, this number is set to 2)
 *
 */
class Interaction
{

  // === PRIVATE MEMBERS ===

private:

  /** name of the Interaction */
  std::string  id;

  /** number specific to each Interaction */
  int number;

  /** size of the interaction, ie size of y[i] and lambda[i] */
  unsigned int interactionSize;

  /** number of relations in the interaction ( equal to interactionSize / nsLawSize ) */
  unsigned int numberOfRelations;

  /** sum of all DS sizes, for DS involved in the interaction */
  unsigned int sizeOfDS;

  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  VectorOfBlocks y;

  /** previous step values for y */
  VectorOfBlocks yOld;

  /** result of the computeInput function */
  VectorOfBlocks lambda;

  /** previous step values for lambda */
  VectorOfBlocks lambdaOld;

  /** the Dynamical Systems concerned by this interaction */
  DynamicalSystemsSet involvedDS;

  /** the Non-smooth Law of the interaction*/
  NonSmoothLaw *nslaw;

  /** the type of Relation of the interaction */
  Relation *relation;

  /** the NonSmoothDynamicalSystem that owns this Interaction */
  NonSmoothDynamicalSystem * NSDS;

  /** the XML object linked to the Interaction to read XML data */
  InteractionXML *interactionxml;

  /** Flags to know if pointers have been allocated inside constructors or not
   *  isXXXAllocatedIn[i][j] = true means that the vector that corresponds to relation j in the derivative i of XXX
   *  has been allocated inside the class.
   */
  AllocationFlags isYAllocatedIn;
  AllocationFlags isYOldAllocatedIn;
  AllocationFlags isLambdaAllocatedIn;
  AllocationFlags isLambdaOldAllocatedIn;
  bool isRelationAllocatedIn;
  bool isNsLawAllocatedIn;

  // === PRIVATE FUNCTIONS ===

  /** default constructor */
  Interaction();

  // === PUBLIC FUNCTIONS ===

public:

  // === CONSTRUCTORS/DESTRUCTOR ===

  /** copy constructor
   *  \param Interaction* : the object to copy
   */
  Interaction(const Interaction& inter);

  /** constructor with XML object of the Interaction
   *  \param InteractionXML* : the XML object corresponding
   *  \param NonSmoothDynamicalSystem (optional)
   */
  Interaction(InteractionXML*, NonSmoothDynamicalSystem* = NULL);

  /** constructor with a set of data
  *  \param string: the id of this Interaction
  *  \param a DynamicalSystemsSet: the set of DS involved in the Interaction
  *  \param int : the number of this Interaction
  *  \param int : the value of interactionSize
  *  \param NonSmoothLaw* : a pointer to the non smooth law
  *  \param Relation* : a pointer to the Relation
  *  \exception RuntimeException
  */
  Interaction(const std::string&, DynamicalSystemsSet&, int, int, NonSmoothLaw*, Relation*);

  /** destructor
   */
  ~Interaction();

  /** allocate memory for y[i] and lambda[i] and set them to zero.
  */
  void initialize();

  /** build Y and Lambda stl vectors.
  *   \param unsigned int: dim of Y and Lambda vector of Interactions (ie number of derivatives that
  *          are taken into account). This is an input argument since it depends on the simulation type.
  */
  void initializeMemory(const unsigned int);

  // === GETTERS/SETTERS ===

  /** get the id of this Interaction
  *  \return the string, id of this Interaction
  */
  inline const std::string  getId() const
  {
    return id;
  }

  /** set the id of this Interaction
   *  \param the integer to set the id
   */
  inline void setId(const int newId)
  {
    id = newId;
  }

  /** get the value of number
   *  \return the value of number
   */
  inline const int getNumber() const
  {
    return number;
  }

  /** set number
  *  \param int number : the value to set number
  */
  inline void setNumber(const int newNumber)
  {
    number = newNumber;
  }

  /** get the dimension of the interaction (y and lambda size)
  *  \return an unsigned int
  */
  inline const unsigned int getInteractionSize() const
  {
    return interactionSize;
  }

  /** set the dimension of the Interaction
  *  \param an unsigned int
  */
  inline void setInteractionSize(const unsigned int newVal)
  {
    interactionSize = newVal;
  }

  /** get the number of relations in the interaction
  *  \return an unsigned int
  */
  inline const unsigned int getNumberOfRelations() const
  {
    return numberOfRelations;
  }

  /** set the number of relations
  *  \param an unsigned int
  */
  inline void setNumberOfRelations(const unsigned int newVal)
  {
    numberOfRelations = newVal;
  }

  /** get the sum of DS sizes, for DS involved in interaction
  *  \return an unsigned int
  */
  inline const unsigned int getSizeOfDS() const
  {
    return sizeOfDS;
  }

  // -- y --

  /** get vector of output derivatives
  *  \return a VectorOfBlocks
  */
  inline const VectorOfBlocks getY() const
  {
    return y;
  }

  /** get y[i], derivative number i of output
  *  \return BlockVector
  */
  inline const BlockVector getY(const unsigned int i) const
  {
    return *(y[i]);
  }

  /** get y[i], derivative number i of output
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getYPtr(const unsigned int i) const
  {
    return y[i];
  }

  /** set the output vector y to newVector with copy of the y[i] (ie memory allocation)
  *  \param VectorOfBlocks
  */
  void setY(const VectorOfBlocks&);

  /** set the output vector y to newVector with direct pointer equality for the y[i]
  *  \param VectorOfBlocks
  */
  void setYPtr(const VectorOfBlocks&);

  /** set y[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setY(const unsigned int , const BlockVector&);

  /** set y[i] to pointer newPtr
  *  \param a SiconosVector * and an unsigned int
  */
  void setYPtr(const unsigned int , SiconosVector *newPtr);

  // -- yOld --

  /** get vector of output derivatives
  *  \return a VectorOfBlocks
  */
  inline const VectorOfBlocks getYOld() const
  {
    return yOld;
  }

  /** get yOld[i], derivative number i of output
  *  \return BlockVector
  */
  inline const BlockVector getYOld(const unsigned int i) const
  {
    return *(yOld[i]);
  }

  /** get yOld[i], derivative number i of output
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getYOldPtr(const unsigned int i) const
  {
    return yOld[i];
  }

  /** set the output vector yOld to newVector
  *  \param VectorOfBlocks
  */
  void setYOld(const VectorOfBlocks&);

  /** set vector yOld to newVector with direct pointer equality for the yOld[i]
  *  \param VectorOfBlocks
  */
  void setYOldPtr(const VectorOfBlocks&);

  /** set yOld[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setYOld(const unsigned int , const BlockVector&);

  /** set yOld[i] to pointer newPtr
  *  \param a SiconosVector * and an unsigned int
  */
  void setYOldPtr(const unsigned int , SiconosVector *newPtr);

  // -- lambda --

  /** get vector of input derivatives
  *  \return a VectorOfBlocks
  */
  inline const VectorOfBlocks getLambda() const
  {
    return lambda;
  }

  /** get lambda[i], derivative number i of input
  *  \return BlockVector
  */
  inline const BlockVector getLambda(const unsigned int i) const
  {
    return *(lambda[i]);
  }

  /** get lambda[i], derivative number i of input
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getLambdaPtr(const unsigned int i) const
  {
    return lambda[i];
  }

  /** set the input vector lambda to newVector
  *  \param VectorOfBlocks
  */
  void setLambda(const VectorOfBlocks&);

  /** set vector lambda to newVector with direct pointer equality for the lambda[i]
  *  \param VectorOfBlocks
  */
  void setLambdaPtr(const VectorOfBlocks&);

  /** set lambda[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setLambda(const unsigned int , const BlockVector&);

  /** set lambda[i] to pointer newPtr
  *  \param a SiconosVector * and an unsigned int
  */
  void setLambdaPtr(const unsigned int , SiconosVector *newPtr);

  // -- lambdaOld --

  /** get vector of input derivatives
  *  \return a VectorOfBlocks
  */
  inline const VectorOfBlocks getLambdaOld() const
  {
    return lambdaOld;
  }

  /** get lambdaOld[i], derivative number i of input
  *  \return BlockVector
  */
  inline const BlockVector getLambdaOld(const unsigned int i) const
  {
    return *(lambdaOld[i]);
  }

  /** get lambdaOld[i], derivative number i of input
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getLambdaOldPtr(const unsigned int i) const
  {
    return lambdaOld[i];
  }

  /** set the input vector lambdaOld to newVector
  *  \param VectorOfBlocks
  */
  void setLambdaOld(const VectorOfBlocks&);

  /** set vector lambdaOld to newVector with direct pointer equality for the lambdaOld[i]
  *  \param VectorOfBlocks
  */
  void setLambdaOldPtr(const VectorOfBlocks&);

  /** set lambdaOld[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setLambdaOld(const unsigned int , const BlockVector&);

  /** set lambdaOld[i] to pointer newPtr
  *  \param a SiconosVector * and an unsigned int
  */
  void setLambdaOldPtr(const unsigned int , SiconosVector *newPtr);

  /** get the DynamicalSystems of this Interaction
  *  \return a DynamicalSystemsSet
  */
  inline DynamicalSystemsSet getDynamicalSystems() const
  {
    return involvedDS;
  }

  /** set the involvedDS
  *  \param a DynamicalSystemsSet
  */
  void setDynamicalSystems(const DynamicalSystemsSet&) ;

  /** get a specific DynamicalSystem
  *  \param the identification number of the wanted DynamicalSystem
  *  \return a pointer on Dynamical System
  */
  DynamicalSystem* getDynamicalSystemPtr(const int);

  /** get a specific DynamicalSystem
  *  \param the identification number of the wanted DynamicalSystem
  *  \return a Dynamical System
  */
  DynamicalSystem getDynamicalSystem(const int);

  /** get the Relation of this Interaction
  *  \return a pointer on this Relation
  */
  inline Relation* getRelationPtr() const
  {
    return relation;
  }

  /** set the Relation of this Interaction
  *  \param the relation* to set
  */
  void setRelationPtr(Relation* newRelation) ;

  /** get the NonSmoothLaw of this Interaction
  *  \return a pointer on this NonSmoothLaw
  */
  inline NonSmoothLaw* getNonSmoothLawPtr() const
  {
    return nslaw;
  }

  /** set the NonSmoothLaw of this Interaction
  *  \param the NonSmoothLaw* to set
  */
  void setNonSmoothLawPtr(NonSmoothLaw* newNslaw) ;

  /** get the NonSmoothDynamicalSystem that contains the current Interaction
  *  \return NonSmoothDynamicalSystem*
  */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return NSDS;
  }

  /** set the NonSmoothDynamicalSystem that contains the current Interaction
  *  \param NonSmoothDynamicalSystem*
  */
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newNsds)
  {
    NSDS = newNsds;
  }

  /** get the InteractionXML* of the Interaction
  *  \return InteractionXML* : the pointer on the InteractionXML
  */
  inline InteractionXML* getInteractionXMLPtr() const
  {
    return interactionxml;
  }

  /** set the InteractionXML* of the Interaction
  *  \param InteractionXML* :  the pointer to set
  */
  inline void setInteractionXMLPtr(InteractionXML* interxml)
  {
    interactionxml = interxml;
  }

  // --- OTHER FUNCTIONS ---

  /** compute sum of all interaction-involved DS sizes
  */
  void computeSizeOfDS();

  /**   put values of y in yOld, the same for lambda
  */
  void swapInMemory();

  /** print the data to the screen
  */
  void display() const;

  // --- XML RELATED FUNCTIONS ---

  /** copy the data of the Interaction to the XML tree
  *  \exception RuntimeException
  */
  void saveInteractionToXML();

};

#endif // INTERACTION_H
