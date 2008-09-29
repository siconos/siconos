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

/*! \file Interaction.h
  \brief Interaction class and related typedef
*/


#ifndef INTERACTION_H
#define INTERACTION_H

// const
#include <boost/shared_ptr.hpp>
#include "BlockVector.h"
#include "DynamicalSystemsSet.h"
#include "Tools.h"

#include "SiconosPointers.h"

class NonSmoothLaw;
class DynamicalSystem;
class Relation;
class NonSmoothDynamicalSystem;
class InteractionXML;
class BlockVector;

/**  An Interaction describes the non-smooth interactions between some Dynamical Systems.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 * An interaction represents the "link" between a set of Dynamical Systems (var: involvedDS) that interact through
 * some relations (between state variables (x,R) and local variables (y,lambda)) completed by a non-smooth law.
 *
 * Thus, the interaction main members are:
 *
 * - a set of Dynamical Systems (from 1 to ...) that interacts, named involvedDS.
 *
 * - relation: a pointer to a Relation object that determines the type of relation and so the way it is computed.\n
 *   Warning: there is only one Relation object (ie only one type of relation for an interaction) but there can be several "relations", in the sense
 *   of constraints equations between (y,lambda) and (x,r).
 *
 * - nslaw: the non smooth law
 *
 * - the local variables y and lambda (their size is interactionSize).
 *   STL vectors are used and y[i] (resp lambda[i]) represents the i-eme derivative of variable y (resp lambda).
 *
 *   y is a container of BlockVector. Each block corresponds to a "unitary relation", a relation which size is the one of the non-smooth law. \n
 *    => ySize = interactionSize = numberOfRelations * nsLawSize .
 *   Same thing for lambda.
 *
 *  => all the relations of the interaction have the same non-smooth law. \n
 *  => the number of relations is equal to interactionSize/nsLawSize. Thus a relation is not necessarily represented
 *     by a single equation.
 *
 *
 */
class Interaction : public boost::enable_shared_from_this<Interaction>
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

  /** sum of all z sizes, for DS involved in the interaction */
  unsigned int sizeZ;

  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  VectorOfVectors y;

  /** previous step values for y */
  VectorOfVectors yOld;

  /** result of the computeInput function */
  VectorOfVectors lambda;

  /** previous step values for lambda */
  VectorOfVectors lambdaOld;

  /** the Dynamical Systems concerned by this interaction */
  DynamicalSystemsSetSPtr involvedDS;

  /** the Non-smooth Law of the interaction*/
  NonSmoothLawSPtr nslaw;

  /** the type of Relation of the interaction */
  RelationSPtr relation;

  /** the NonSmoothDynamicalSystem that owns this Interaction */
  SP::NonSmoothDynamicalSystem NSDS;

  /** the XML object linked to the Interaction to read XML data */
  InteractionXMLSPtr interactionxml;

  // === PRIVATE FUNCTIONS ===

  /** default constructor */
  Interaction();

  /** copy constructor => private, no copy nor pass-by-value.
   */
  Interaction(const Interaction& inter);

public:

  /** constructor with XML object of the Interaction
   *  \param InteractionXML* : the XML object corresponding
   *  \param NonSmoothDynamicalSystem (optional)
   */
  Interaction(InteractionXMLSPtr, SP::NonSmoothDynamicalSystem = boost::shared_ptr<NonSmoothDynamicalSystem>());

  /** constructor with a set of data (only one DS in the Interaction) - Note: no id.
   *  \param a SP::DynamicalSystem: the DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of interactionSize
   *  \param NonSmoothLaw* : a pointer to the non smooth law
   *  \param Relation* : a pointer to the Relation
   */
  Interaction(SP::DynamicalSystem, int, int, NonSmoothLawSPtr, RelationSPtr);
  /** constructor with a set of data (only one DS in the Interaction)
   *  \param string: the id of this Interaction
   *  \param a SP::DynamicalSystem: the DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of interactionSize
   *  \param NonSmoothLaw* : a pointer to the non smooth law
   *  \param Relation* : a pointer to the Relation
   */
  Interaction(const std::string&, SP::DynamicalSystem, int, int, NonSmoothLawSPtr, RelationSPtr);

  /** constructor with a set of data - Note: no id.
   *  \param a DynamicalSystemsSet: the set of DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of interactionSize
   *  \param NonSmoothLaw* : a pointer to the non smooth law
   *  \param Relation* : a pointer to the Relation
   */
  Interaction(DynamicalSystemsSet&, int, int, NonSmoothLawSPtr, RelationSPtr);

  /** constructor with a set of data
   *  \param string: the id of this Interaction
   *  \param a DynamicalSystemsSet: the set of DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of interactionSize
   *  \param NonSmoothLaw* : a pointer to the non smooth law
   *  \param Relation* : a pointer to the Relation
   */
  Interaction(const std::string&, DynamicalSystemsSet&, int, int, NonSmoothLawSPtr, RelationSPtr);

  /** destructor
   */
  ~Interaction();

  /** allocate memory for y[i] and lambda[i] and set them to zero.
   * \param time for initialization.
   * \param the number of required derivatives for y.
   */
  void initialize(double, unsigned int);

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
  inline const unsigned int getSizeOfY() const
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

  /** get the sum of DS sizes, for DS involved in interaction
   *  \return an unsigned int
   */
  inline const unsigned int getSizeOfDS() const
  {
    return sizeOfDS;
  }

  /** get the sum of z sizes, for DS involved in interaction
   *  \return an unsigned int
   */
  inline const unsigned int getSizeZ() const
  {
    return sizeZ;
  }

  // -- y --

  /** get vector of output derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getY() const
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
  inline SiconosVectorSPtr getYPtr(const unsigned int i) const
  {
    return y[i];
  }

  /** set the output vector y to newVector with copy of the y[i] (ie memory allocation)
  *  \param VectorOfVectors
  */
  void setY(const VectorOfVectors&);

  /** set the output vector y to newVector with direct pointer equality for the y[i]
  *  \param VectorOfVectors
  */
  void setYPtr(const VectorOfVectors&);

  /** set y[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setY(const unsigned int , const BlockVector&);

  /** set y[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setYPtr(const unsigned int , SiconosVectorSPtr newPtr);

  // -- yOld --

  /** get vector of output derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getYOld() const
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
  inline SiconosVectorSPtr getYOldPtr(const unsigned int i) const
  {
    return yOld[i];
  }

  /** set the output vector yOld to newVector
  *  \param VectorOfVectors
  */
  void setYOld(const VectorOfVectors&);

  /** set vector yOld to newVector with direct pointer equality for the yOld[i]
  *  \param VectorOfVectors
  */
  void setYOldPtr(const VectorOfVectors&);

  /** set yOld[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setYOld(const unsigned int , const BlockVector&);

  /** set yOld[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setYOldPtr(const unsigned int , SiconosVectorSPtr newPtr);

  // -- lambda --

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getLambda() const
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
  inline SiconosVectorSPtr getLambdaPtr(const unsigned int i) const
  {
    return lambda[i];
  }

  /** set the input vector lambda to newVector
  *  \param VectorOfVectors
  */
  void setLambda(const VectorOfVectors&);

  /** set vector lambda to newVector with direct pointer equality for the lambda[i]
  *  \param VectorOfVectors
  */
  void setLambdaPtr(const VectorOfVectors&);

  /** set lambda[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setLambda(const unsigned int , const BlockVector&);

  /** set lambda[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setLambdaPtr(const unsigned int , SiconosVectorSPtr newPtr);

  // -- lambdaOld --

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getLambdaOld() const
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
  inline SiconosVectorSPtr getLambdaOldPtr(const unsigned int i) const
  {
    return lambdaOld[i];
  }

  /** set the input vector lambdaOld to newVector
  *  \param VectorOfVectors
  */
  void setLambdaOld(const VectorOfVectors&);

  /** set vector lambdaOld to newVector with direct pointer equality for the lambdaOld[i]
  *  \param VectorOfVectors
  */
  void setLambdaOldPtr(const VectorOfVectors&);

  /** set lambdaOld[i] to newValue
  *  \param a BlockVector and an unsigned int
  */
  void setLambdaOld(const unsigned int , const BlockVector&);

  /** set lambdaOld[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setLambdaOldPtr(const unsigned int , SiconosVectorSPtr newPtr);

  /** gets an iterator to the first element of the involvedDS set.
   *  \return a DSIterator.
   */

  inline DSIterator dynamicalSystemsBegin()
  {
    return involvedDS->begin();
  };

  /** gets an iterator equal to involvedDS->end().
   *  \return a DSIterator.
   */
  inline DSIterator dynamicalSystemsEnd()
  {
    return involvedDS->end();
  };

  /** gets a const iterator to the first element of the involvedDS set.
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsBegin() const
  {
    return involvedDS->begin();
  };

  /** gets a const iterator equal to involvedDS->end().
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsEnd() const
  {
    return involvedDS->end();
  };

  /** get a pointer to the DynamicalSystems of this Interaction
   *  \return a DynamicalSystemsSet*
   */
  inline DynamicalSystemsSetSPtr getDynamicalSystemsPtr()
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
  SP::DynamicalSystem getDynamicalSystemPtr(int);

  /** get a specific DynamicalSystem. Out of date function?
   *  \param the identification number of the wanted DynamicalSystem
   *  \param a Dynamical System: out-parameter
   */
  void getDynamicalSystem(int, DynamicalSystem&);

  /** get the Relation of this Interaction
  *  \return a pointer on this Relation
  */
  inline RelationSPtr getRelationPtr() const
  {
    return relation;
  }

  /** set the Relation of this Interaction
  *  \param the relation* to set
  */
  void setRelationPtr(RelationSPtr newRelation) ;

  /** get the NonSmoothLaw of this Interaction
  *  \return a pointer on this NonSmoothLaw
  */
  inline NonSmoothLawSPtr getNonSmoothLawPtr() const
  {
    return nslaw;
  }

  /** set the NonSmoothLaw of this Interaction
  *  \param the NonSmoothLaw* to set
  */
  void setNonSmoothLawPtr(NonSmoothLawSPtr newNslaw) ;

  /** get the NonSmoothDynamicalSystem that contains the current Interaction
  *  \return SP::NonSmoothDynamicalSystem
  */
  inline SP::NonSmoothDynamicalSystem getNonSmoothDynamicalSystemPtr() const
  {
    return NSDS;
  }

  /** set the NonSmoothDynamicalSystem that contains the current Interaction
  *  \param SP::NonSmoothDynamicalSystem
  */
  inline void setNonSmoothDynamicalSystemPtr(boost::shared_ptr<NonSmoothDynamicalSystem> newNsds)
  {
    NSDS = newNsds;
  }

  /** function used to sort Interaction in SiconosSet<SP::Interaction>
   *  \return an int
   */
  inline double* const getSort() const
  {
    return (double*)this;
  }

  // --- OTHER FUNCTIONS ---

  /** compute sum of all interaction-involved DS sizes
  */
  void computeSizeOfDS();

  /**   put values of y into yOld, the same for lambda
  */
  void swapInMemory();

  /** print the data to the screen
  */
  void display() const;

  /** Computes output y; depends on the relation type.
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** Compute input r of all Dynamical Systems involved in the present Interaction.
   *  \param double : current time
   *  \param unsigned int: order of lambda used to compute input.
   */
  void computeInput(double, unsigned int);

  // --- XML RELATED FUNCTIONS ---

  /** get the InteractionXML* of the Interaction
  *  \return InteractionXML* : the pointer on the InteractionXML
  */
  inline InteractionXMLSPtr getInteractionXMLPtr() const
  {
    return interactionxml;
  }

  /** set the InteractionXML* of the Interaction
  *  \param InteractionXML* :  the pointer to set
  */
  inline void setInteractionXMLPtr(boost::shared_ptr<InteractionXML> interxml)
  {
    interactionxml = interxml;
  }

  /** copy the data of the Interaction to the XML tree
  *  \exception RuntimeException
  */
  void saveInteractionToXML();

};

#endif // INTERACTION_H
