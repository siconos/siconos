#ifndef INTERACTION_H
#define INTERACTION_H

// for composition ...
#include "DynamicalSystem.h"
#include "NonSmoothLaw.h"
#include "Relation.h"
#include "RelationXML.h"
#include "NonSmoothDynamicalSystem.h"

// IO (XML)
#include "InteractionXML.h"

// tools
#include "SimpleVector.h"

// const
#include "SiconosConst.h"

#include <vector>
#include <string>

class NonSmoothLaw;
class DynamicalSystem;
class Relation;
class NonSmoothDynamicalSystem;

/** \class Interaction
 *  \brief this class describes interaction between some dynamical systems (DS) (from 1 to NumberOfDS)
 *  Each interaction owns local variables  with relations between them, plus a non smooth law
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 * \todo save a number of derivatives of y depending on the relative degree
 * \todo remove createRelay etc ...
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

  /** number of relations in the interaction (ie size of y[i], lambda[i] ...*/
  unsigned int nInteraction;

  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  std::vector< SimpleVector* >  y;

  /** previous step values for y */
  std::vector< SimpleVector* >  yOld;

  /** result of the computeInput function */
  std::vector< SimpleVector* >  lambda;

  /** previous step values for lambda */
  std::vector< SimpleVector* > lambdaOld;

  /** the Dynamical Systems concerned by this interaction
   *  their number is between 1 and NumberOfDs */
  std::vector<DynamicalSystem*> vectorDS;

  /** the Non-smooth Law of the interaction*/
  NonSmoothLaw *nslaw;

  /** the type of Relation of the interaction */
  Relation *relation;

  /** the XML object linked to the Interaction to read XML data */
  InteractionXML *interactionxml;

  /** Flags to know if pointers have been allocated inside constructors or not */
  std::vector<bool> isYAllocatedIn;
  std::vector<bool> isYOldAllocatedIn;
  std::vector<bool> isLambdaAllocatedIn;
  std::vector<bool> isLambdaOldAllocatedIn;
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

  /** \fn Interaction(const string&, const int& number,const int& nInter,
   *    vector<DynamicalSystem*>* dsConcerned)
   *  \brief constructor with a set of data
   *  \param string: the id of this Interaction
   *  \param int : the number of this Interaction (optional)
   *  \param int : the value of nInter for this Interaction (optional)
   *  \param vector<DynamicalSystem*>* : the Dynamical Systems concerned by this interaction (optional)
   *  \exception RuntimeException
   */
  Interaction(const std::string&, const int& = -1 , const int& = -1, std::vector<DynamicalSystem*>* = NULL);

  /** \fn ~Interaction()
   * \brief destructor
   */
  ~Interaction();

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

  /** \fn void setNumber(const int&)
   *  \brief set number
   *  \param int number : the value to set number
   */
  inline void setNumber(const int& newNumber)
  {
    number = newNumber;
  }

  /** \fn const int getNInteraction() const
   *  \brief get the number of relations of this Interaction
   *  \return the value of nInteraction
   */
  inline const unsigned int getNInteraction() const
  {
    return nInteraction;
  }

  /** \fn void setNInteraction(const int&)
   *  \brief set the dimension of the Interaction
   *  \param the integer to set the nInteraction
   */
  inline void setNInteraction(const unsigned int& nInter)
  {
    nInteraction = nInter;
  }

  // -- y --

  /** \fn  const vector<SimpleVector*> getY() const
   *  \brief get vector of output derivatives
   *  \return a vector<SimpleVector*>
   */
  inline const std::vector<SimpleVector*> getY() const
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

  /** \fn void setY (const vector<SimpleVector*>& newVector)
   *  \brief set the output vector y to newVector with copy of the y[i] (ie memory allocation)
   *  \param std::vector<SimpleVector*>
   */
  void setY(const std::vector<SimpleVector*>&);

  /** \fn void setYPtr (const vector<SimpleVector*>& newVector)
   *  \brief set the output vector y to newVector with direct pointer equality for the y[i]
   *  \param std::vector<SimpleVector*>
   */
  void setYPtr(const std::vector<SimpleVector*>&);

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

  /** \fn  const vector<SimpleVector*> getYOld() const
   *  \brief get vector of output derivatives
   *  \return a vector<SimpleVector*>
   */
  inline const std::vector<SimpleVector*> getYOld() const
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

  /** \fn void setYOld (const vector<SimpleVector*>& newVector)
   *  \brief set the output vector yOld to newVector
   *  \param std::vector<SimpleVector*>
   */
  void setYOld(const std::vector<SimpleVector*>&);

  /** \fn void setYOldPtr(const std::vector<SimpleVector*>& newVector);
   *  \brief set vector yOld to newVector with direct pointer equality for the yOld[i]
   *  \param std::vector<SimpleVector*>
   */
  void setYOldPtr(const std::vector<SimpleVector*>&);

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

  /** \fn  const vector<SimpleVector*> getLambda() const
   *  \brief get vector of input derivatives
   *  \return a vector<SimpleVector*>
   */
  inline const std::vector<SimpleVector*> getLambda() const
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

  /** \fn void setLambda (const vector<SimpleVector*>& newVector)
   *  \brief set the input vector lambda to newVector
   *  \param std::vector<SimpleVector*>
   */
  void setLambda(const std::vector<SimpleVector*>&);

  /** \fn void setLambdaPtr(const std::vector<SimpleVector*>& newVector);
   *  \brief set vector lambda to newVector with direct pointer equality for the lambda[i]
   *  \param std::vector<SimpleVector*>
   */
  void setLambdaPtr(const std::vector<SimpleVector*>&);

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

  /** \fn  const vector<SimpleVector*> getLambdaOld() const
   *  \brief get vector of input derivatives
   *  \return a vector<SimpleVector*>
   */
  inline const std::vector<SimpleVector*> getLambdaOld() const
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

  /** \fn void setLambdaOld (const vector<SimpleVector*>& newVector)
   *  \brief set the input vector lambdaOld to newVector
   *  \param std::vector<SimpleVector*>
   */
  void setLambdaOld(const std::vector<SimpleVector*>&);

  /** \fn void setLambdaOldPtr(const std::vector<SimpleVector*>& newVector);
   *  \brief set vector lambdaOld to newVector with direct pointer equality for the lambdaOld[i]
   *  \param std::vector<SimpleVector*>
   */
  void setLambdaOldPtr(const std::vector<SimpleVector*>&);

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

  /** \fn vector<DynamicalSystem*> getDynamicalSystems()
   *  \brief get the DynamicalSystems of this Interaction
   *  \return vector<DynamicalSystem*> : a vector of DS
   */
  inline std::vector<DynamicalSystem*> getDynamicalSystems() const
  {
    return vectorDS;
  }

  /** \fn void setDynamicalSystems(const std::vector<DynamicalSystem*>&)
   *  \brief set the ds vector
   *  \param a std vector of DynamicalSystem*>
   */
  void setDynamicalSystems(const std::vector<DynamicalSystem*>&) ;

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

  /** \fn NonSmoothLaw* createNewtonImpactLawNSL(const double& e)
   *  \brief allows to create a NewtonImpactLawNSL non-smooth law for this Interaction
   *  \param double : the e value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createNewtonImpactLawNSL(const double&);

  /** \fn NonSmoothLaw* createNewtonImpactLawNSL(const double& en, const double& et, const double& mu)
   *  \brief allows to create a NewtonImpactLawNSL non-smooth law for this Interaction
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
