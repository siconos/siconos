/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#ifndef INTERACTION_H
#define INTERACTION_H

// for composition ...
#include "DynamicalSystem.h"
#include "NonSmoothLaw.h"
#include "Relation.h"
#include "RelationXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "DSSet.h"

// IO (XML)
#include "InteractionXML.h"

// tools
#include "SimpleVector.h"

// const
#include "SiconosConst.h"

// stl tools
#include <vector>
#include <string>
#include <set>

class DSSet;
class NonSmoothLaw;
class DynamicalSystem;
class Relation;
class NonSmoothDynamicalSystem;

/** \class Interaction
 *  \brief this class describes interaction between some dynamical systems (DS) (from 1 to NumberOfDS)
 *  Each interaction owns local variables  with relations between them, plus a non smooth law
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) Apr 29, 2004
 *
 * An interaction represents the "link" between a set of Dynamical Systems (var: involvedDS) that interact alltogether through
 * some relations (between state variables (x,R) and local ones (y,lambda)) completed by a non-smooth law.
 *
 * Thus, the interaction main members are:
 * - a set of Dynamical Systems that interacts, involvedDS.
 * - y and lambda (their size is interactionSize).
 *   stl vectors are used and y[i] (resp lambda[i]) represents the i-eme derivative of variable y (resp lambda)
 * - relation: a pointer to a Relation object that determines the type of relation and so the way it is computed.
 * - nslaw: the non smooth law
 *
 * Note that numberOfRelations is equal to interactionSize/nsLawSize. Thus a relation is not necessarily represented
 * by a single equation.
 *
 * \todo save a number of derivatives of y depending on the relative degree (at the time, this number is set to 2)
 * \todo remove createRelay etc ...
 *
 */

/** container for SiconosVectors*/
typedef std::vector< SimpleVector* > vectorOfSiconosVector ;

/** iterator through vector of SimpleVector */
typedef vectorOfSiconosVector::iterator vectorOfSiconosVectorIt;

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
  vectorOfSiconosVector y;

  /** previous step values for y */
  vectorOfSiconosVector yOld;

  /** result of the computeInput function */
  vectorOfSiconosVector lambda;

  /** previous step values for lambda */
  vectorOfSiconosVector lambdaOld;

  /** the Dynamical Systems concerned by this interaction */
  DSSet involvedDS;

  /** the Non-smooth Law of the interaction*/
  NonSmoothLaw *nslaw;

  /** the type of Relation of the interaction */
  Relation *relation;

  /** the NonSmoothDynamicalSystem that owns this Interaction */
  NonSmoothDynamicalSystem * NSDS;

  /** the XML object linked to the Interaction to read XML data */
  InteractionXML *interactionxml;

  /** Flags to know if pointers have been allocated inside constructors or not */
  std::deque<bool> isYAllocatedIn;
  std::deque<bool> isYOldAllocatedIn;
  std::deque<bool> isLambdaAllocatedIn;
  std::deque<bool> isLambdaOldAllocatedIn;
  bool isRelationAllocatedIn;
  bool isNsLawAllocatedIn;

  // === PRIVATE FUNCTIONS ===

  /** \fn Interaction()
   *  \brief default constructor
   */
  Interaction();

  // === PUBLIC FUNCTIONS ===

public:

  // === CONSTRUCTORS/DESTRUCTOR ===

  /** \fn Interaction(const Interaction&)
   *  \brief copy constructor
   *  \param Interaction* : the object to copy
   */
  Interaction(const Interaction& inter);

  /** \fn Interaction(InteractionXML*)
   *  \brief constructor with XML object of the Interaction
   *  \param InteractionXML* : the XML object corresponding
   *  \param NonSmoothDynamicalSystem : the nsds that owns this strategy (optional)
   */
  Interaction(InteractionXML*, NonSmoothDynamicalSystem* = NULL);

  /** \fn Interaction(const string&, DSSet& dsConcerned, const int& number,const int& nInter)
   *  \brief constructor with a set of data
   *  \param a DSSet: the set of DS involved in the Interaction
   *  \param string: the id of this Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of nInter for this Interaction (optional)
   *  \exception RuntimeException
   */
  Interaction(const std::string&, DSSet&, const int&, const int& = -1);

  /** \fn ~Interaction()
   * \brief destructor
   */
  ~Interaction();

  /** \fn void initialize()
   *  \brief allocate memory for y[i] and lambda[i] and set them to zero.
   */
  void initialize();

  /** \fn void initializeVectors(vectorOfSiconosVector inputVector)
   *  \brief set all components of inputVector to zero
   *  \param a vector of SimpleVector*
   */
  void initializeVectors(vectorOfSiconosVector);

  // === GETTERS/SETTERS ===

  /** \fn const string getId() const
   *  \brief get the id of this Interaction
   *  \return the string, id of this Interaction
   */
  inline const std::string  getId() const
  {
    return id;
  }

  /** \fn void setId(int)
   *  \brief set the id of this Interaction
   *  \param the integer to set the id
   */
  inline void setId(const int& newId)
  {
    id = newId;
  }

  /** \fn const int getNumber() const
   *  \brief get the value of number
   *  \return the value of number
   */
  inline const int getNumber() const
  {
    return number;
  }

  /** \fn const int getNumberForSorting(void) const;
   *  \brief same as getNumber, but return an unsigned long int, used for set<Interactions*> in OSI, NSDS ...
   *   as sorting criterion.
   *  \return the value of number
   */
  inline const unsigned long int getNumberForSorting() const
  {
    return number;
  }

  /** \fn void setNumber(const int&)
   *  \brief set number
   *  \param int number : the value to set number
   */
  inline void setNumber(const int& newNumber)
  {
    number = newNumber;
  }

  /** \fn const unsigned int getInteractionSize() const
   *  \brief get the dimension of the interaction (y and lambda size)
   *  \return an unsigned int
   */
  inline const unsigned int getInteractionSize() const
  {
    return interactionSize;
  }

  /** \fn void setInteractionSize(const unsigned int&)
   *  \brief set the dimension of the Interaction
   *  \param an unsigned int
   */
  inline void setInteractionSize(const unsigned int& newVal)
  {
    interactionSize = newVal;
  }

  /** \fn const unsigned int getNumberOfRelations() const
   *  \brief get the number of relations in the interaction
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfRelations() const
  {
    return numberOfRelations;
  }

  /** \fn void setNumberOfRelations(const unsigned int&)
   *  \brief set the number of relations
   *  \param an unsigned int
   */
  inline void setNumberOfRelations(const unsigned int& newVal)
  {
    numberOfRelations = newVal;
  }

  /** \fn const unsigned int getSizeOfDS() const
   *  \brief get the sum of DS sizes, for DS involved in interaction
   *  \return an unsigned int
   */
  inline const unsigned int getSizeOfDS() const
  {
    return sizeOfDS;
  }

  // -- y --

  /** \fn  const vectorOfSiconosVector getY() const
   *  \brief get vector of output derivatives
   *  \return a vectorOfSiconosVector
   */
  inline const vectorOfSiconosVector getY() const
  {
    return y;
  }

  /** \fn  const SimpleVector getY(const unsigned int & i) const
   *  \brief get y[i], derivative number i of output
   *  \return SimpleVector
   */
  inline const SimpleVector getY(const unsigned int& i) const
  {
    return *(y[i]);
  }

  /** \fn SimpleVector* getYPtr(const unsigned int& i) const
   *  \brief get y[i], derivative number i of output
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYPtr(const unsigned int& i) const
  {
    return y[i];
  }

  /** \fn void setY (const vectorOfSiconosVector& newVector)
   *  \brief set the output vector y to newVector with copy of the y[i] (ie memory allocation)
   *  \param vectorOfSiconosVector
   */
  void setY(const vectorOfSiconosVector&);

  /** \fn void setYPtr (const vectorOfSiconosVector& newVector)
   *  \brief set the output vector y to newVector with direct pointer equality for the y[i]
   *  \param vectorOfSiconosVector
   */
  void setYPtr(const vectorOfSiconosVector&);

  /** \fn void setY (const unsigned int & i, const SimpleVector& newValue);
   *  \brief set y[i] to newValue
   *  \param a SimpleVector and an unsigned int
   */
  void setY(const unsigned int &, const SimpleVector&);

  /** \fn void setYPtr(const unsigned int & i, SimpleVector* newPtr)
   *  \brief set y[i] to pointer newPtr
   *  \param a SimpleVector * and an unsigned int
   */
  void setYPtr(const unsigned int &, SimpleVector *newPtr);

  // -- yOld --

  /** \fn  const vectorOfSiconosVector getYOld() const
   *  \brief get vector of output derivatives
   *  \return a vectorOfSiconosVector
   */
  inline const vectorOfSiconosVector getYOld() const
  {
    return yOld;
  }

  /** \fn  const SimpleVector getYOld(const unsigned int & i) const
   *  \brief get yOld[i], derivative number i of output
   *  \return SimpleVector
   */
  inline const SimpleVector getYOld(const unsigned int& i) const
  {
    return *(yOld[i]);
  }

  /** \fn SimpleVector* getYOldPtr(const unsigned int& i) const
   *  \brief get yOld[i], derivative number i of output
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYOldPtr(const unsigned int& i) const
  {
    return yOld[i];
  }

  /** \fn void setYOld (const vectorOfSiconosVector& newVector)
   *  \brief set the output vector yOld to newVector
   *  \param vectorOfSiconosVector
   */
  void setYOld(const vectorOfSiconosVector&);

  /** \fn void setYOldPtr(const vectorOfSiconosVector& newVector);
   *  \brief set vector yOld to newVector with direct pointer equality for the yOld[i]
   *  \param vectorOfSiconosVector
   */
  void setYOldPtr(const vectorOfSiconosVector&);

  /** \fn void setYOld (const unsigned int & i, const SimpleVector& newValue);
   *  \brief set yOld[i] to newValue
   *  \param a SimpleVector and an unsigned int
   */
  void setYOld(const unsigned int &, const SimpleVector&);

  /** \fn void setYOldPtr(const unsigned int & i, SimpleVector* newPtr)
   *  \brief set yOld[i] to pointer newPtr
   *  \param a SimpleVector * and an unsigned int
   */
  void setYOldPtr(const unsigned int &, SimpleVector *newPtr);

  // -- lambda --

  /** \fn  const vectorOfSiconosVector getLambda() const
   *  \brief get vector of input derivatives
   *  \return a vectorOfSiconosVector
   */
  inline const vectorOfSiconosVector getLambda() const
  {
    return lambda;
  }

  /** \fn  const SimpleVector getLambda(const unsigned int & i) const
   *  \brief get lambda[i], derivative number i of input
   *  \return SimpleVector
   */
  inline const SimpleVector getLambda(const unsigned int& i) const
  {
    return *(lambda[i]);
  }

  /** \fn SimpleVector* getLambdaPtr(const unsigned int& i) const
   *  \brief get lambda[i], derivative number i of input
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLambdaPtr(const unsigned int& i) const
  {
    return lambda[i];
  }

  /** \fn void setLambda (const vectorOfSiconosVector& newVector)
   *  \brief set the input vector lambda to newVector
   *  \param vectorOfSiconosVector
   */
  void setLambda(const vectorOfSiconosVector&);

  /** \fn void setLambdaPtr(const vectorOfSiconosVector& newVector);
   *  \brief set vector lambda to newVector with direct pointer equality for the lambda[i]
   *  \param vectorOfSiconosVector
   */
  void setLambdaPtr(const vectorOfSiconosVector&);

  /** \fn void setLambda (const unsigned int & i, const SimpleVector& newValue);
   *  \brief set lambda[i] to newValue
   *  \param a SimpleVector and an unsigned int
   */
  void setLambda(const unsigned int &, const SimpleVector&);

  /** \fn void setLambdaPtr(const unsigned int & i, SimpleVector* newPtr)
   *  \brief set lambda[i] to pointer newPtr
   *  \param a SimpleVector * and an unsigned int
   */
  void setLambdaPtr(const unsigned int &, SimpleVector *newPtr);

  // -- lambdaOld --

  /** \fn  const vectorOfSiconosVector getLambdaOld() const
   *  \brief get vector of input derivatives
   *  \return a vectorOfSiconosVector
   */
  inline const vectorOfSiconosVector getLambdaOld() const
  {
    return lambdaOld;
  }

  /** \fn  const SimpleVector getLambdaOld(const unsigned int & i) const
   *  \brief get lambdaOld[i], derivative number i of input
   *  \return SimpleVector
   */
  inline const SimpleVector getLambdaOld(const unsigned int& i) const
  {
    return *(lambdaOld[i]);
  }

  /** \fn SimpleVector* getLambdaOldPtr(const unsigned int& i) const
   *  \brief get lambdaOld[i], derivative number i of input
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLambdaOldPtr(const unsigned int& i) const
  {
    return lambdaOld[i];
  }

  /** \fn void setLambdaOld (const vectorOfSiconosVector& newVector)
   *  \brief set the input vector lambdaOld to newVector
   *  \param vectorOfSiconosVector
   */
  void setLambdaOld(const vectorOfSiconosVector&);

  /** \fn void setLambdaOldPtr(const vectorOfSiconosVector& newVector);
   *  \brief set vector lambdaOld to newVector with direct pointer equality for the lambdaOld[i]
   *  \param vectorOfSiconosVector
   */
  void setLambdaOldPtr(const vectorOfSiconosVector&);

  /** \fn void setLambdaOld (const unsigned int & i, const SimpleVector& newValue);
   *  \brief set lambdaOld[i] to newValue
   *  \param a SimpleVector and an unsigned int
   */
  void setLambdaOld(const unsigned int &, const SimpleVector&);

  /** \fn void setLambdaOldPtr(const unsigned int & i, SimpleVector* newPtr)
   *  \brief set lambdaOld[i] to pointer newPtr
   *  \param a SimpleVector * and an unsigned int
   */
  void setLambdaOldPtr(const unsigned int &, SimpleVector *newPtr);

  /** \fn DSSet getDynamicalSystems()
   *  \brief get the DynamicalSystems of this Interaction
   *  \return a DSSet
   */
  inline DSSet getDynamicalSystems() const
  {
    return involvedDS;
  }

  /** \fn void setDynamicalSystems(const DSSet&)
   *  \brief set the involvedDS
   *  \param a DSSet
   */
  void setDynamicalSystems(const DSSet&) ;

  /** \fn DynamicalSystem* getDynamicalSystemPtr(const int&)
   *  \brief get a specific DynamicalSystem
   *  \param the identification number of the wanted DynamicalSystem
   *  \return a pointer on Dynamical System
   */
  DynamicalSystem* getDynamicalSystemPtr(const int&);

  /** \fn DynamicalSystem getDynamicalSystemNumber(const int& nb)
   *  \brief get a specific DynamicalSystem
   *  \param the identification number of the wanted DynamicalSystem
   *  \return a Dynamical System
  */
  DynamicalSystem getDynamicalSystem(const int&);

  /** \fn Relation* getRelationPtr(void)
   *  \brief get the Relation of this Interaction
   *  \return a pointer on this Relation
   */
  inline Relation* getRelationPtr() const
  {
    return relation;
  }

  /** \fn void setRelationPtr(Relation*)
   *  \brief set the Relation of this Interaction
   *  \param the relation* to set
   */
  void setRelationPtr(Relation* newRelation) ;

  /** \fn NonSmoothLaw* getNonSmoothLawPtr(void)
   *  \brief get the NonSmoothLaw of this Interaction
   *  \return a pointer on this NonSmoothLaw
   */
  inline NonSmoothLaw* getNonSmoothLawPtr() const
  {
    return nslaw;
  }

  /** \fn void setNonSmoothLawPtr(NonSmoothLaw*)
   *  \brief set the NonSmoothLaw of this Interaction
   *  \param the NonSmoothLaw* to set
   */
  void setNonSmoothLawPtr(NonSmoothLaw* newNslaw) ;

  /** \fn NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const;
   *  \brief get the NonSmoothDynamicalSystem that contains the current Interaction
   *  \return NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return NSDS;
  }

  /** \fn void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem*);
   *  \brief set the NonSmoothDynamicalSystem that contains the current Interaction
   *  \param NonSmoothDynamicalSystem*
   */
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newNsds)
  {
    NSDS = newNsds;
  }

  /** \fn inline InteractionXML* getInteractionXMLPtr()
   *  \brief get the InteractionXML* of the Interaction
   *  \return InteractionXML* : the pointer on the InteractionXML
   */
  inline InteractionXML* getInteractionXMLPtr() const
  {
    return interactionxml;
  }

  /** \fn inline void setInteractionXMLPtr(InteractionXML* interxml)
   *  \brief set the InteractionXML* of the Interaction
   *  \param InteractionXML* :  the pointer to set
   */
  inline void setInteractionXMLPtr(InteractionXML* interxml)
  {
    interactionxml = interxml;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void computeSizeOfDS()
   * \brief compute sum of all interaction-involved DS sizes
   */
  void computeSizeOfDS();

  /** \fn   void swapInMemory(void);
   * \brief   put values of y in yOld, the same for lambda
   */
  void swapInMemory();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn NonSmoothLaw* createComplementarityConditionNSL()
   *  \brief allows to create a ComplementarityConditionNSL non-smooth law for this Interaction
   *  \return NonSmoothLaw* : the non-smooth law allocated
   */
  NonSmoothLaw* createComplementarityConditionNSL();

  /** \fn NonSmoothLaw* createRelayNSL(const double& c, const double& d)
   *  \brief allows to create a RelayNSL non-smooth law for this Interaction
   *  \param double : the c value
   *  \param double : the d value
   *  \return NonSmoothLaw* : the non-smooth law allocated
   */
  NonSmoothLaw* createRelayNSL(const double&, const double&);

  /** \fn NonSmoothLaw* createNewtonImpactNSL(const double& e)
   *  \brief allows to create a NewtonImpactNSL non-smooth law for this Interaction
   *  \param double : the e value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createNewtonImpactNSL(const double&);

  /** \fn NonSmoothLaw* createNewtonImpactNSL(const double& en, const double& et, const double& mu)
   *  \brief allows to create a NewtonImpactNSL non-smooth law for this Interaction
   *  \param double : the en value
   *  \param double : the et value
   *  \param double : the mu value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createNewtonImpactFrictionNSL(const double&, const double&, const double&);

  // --- XML RELATED FUNCTIONS ---

  /** \fn void saveInteractionToXML()
   *  \brief copy the data of the Interaction to the XML tree
   *  \exception RuntimeException
   */
  void saveInteractionToXML();

};

#endif // INTERACTION_H
