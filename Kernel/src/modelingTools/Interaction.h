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
#include "NewSiconosVector.h"
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
 *  Each interaction owns local variables and some relations between them, plus a non smooth law
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 */
class Interaction
{
public:

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

  /** \fn Interaction(const string&, const int& number,const int& nInter, vector<int>* status,
      vector<DynamicalSystem*>* dsConcerned)
  *  \brief constructor with a set of data
  *  \param string: the id of this Interaction
  *  \param int : the number of this Interaction (optional)
  *  \param int : the value of nInter for this Interaction (optional)
  *  \param vector<int>* : the status of this Interaction (optional)
  *  \param vector<DynamicalSystem*>* : the Dynamical Systems concerned by this interaction (optional)
  *  \exception RuntimeException
  */
  Interaction(const std::string&, const int& = -1 , const int& = -1, std::vector<int>* = NULL, std::vector<DynamicalSystem*>* = NULL);

  ~Interaction();

  // --- GETTERS/SETTERS ---

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
   *  \brief get the number nInteraction of this Interaction
   *  \return the value of nInteraction
   */
  inline const int getNInteraction() const
  {
    return nInteraction;
  }

  /** \fn void setNInteraction(const int&)
   *  \brief allows to set the value of nInteraction
   *  \param the integer to set the nInteraction
   */
  inline void setNInteraction(const int& nInter)
  {
    nInteraction = nInter;
  }

  // -- y --

  /** \fn  const SimpleVector getY() const
   *  \brief get the value of y
   *  \return SimpleVector
   */
  inline const SimpleVector getY() const
  {
    return *y;
  }

  /** \fn SimpleVector* getYPtr() const
   *  \brief get y
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYPtr() const
  {
    return y;
  }

  /** \fn void setY (const SimpleVector& newValue)
   *  \brief set the value of y to newValue
   *  \param SimpleVector newValue
   */
  inline void setY(const SimpleVector& newValue)
  {
    *y = newValue;
  }

  /** \fn void setYPtr(SimpleVector* newPtr)
   *  \brief set Y to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setYPtr(SimpleVector *newPtr)
  {
    if (isYAllocatedIn) delete y;
    y = newPtr;
    isYAllocatedIn = false;
  }

  // -- yDot --

  /** \fn  const SimpleVector getYDot() const
   *  \brief get the value of yDot
   *  \return SimpleVector
   */
  inline const SimpleVector getYDot() const
  {
    return *yDot;
  }

  /** \fn SimpleVector* getYDotPtr() const
   *  \brief get yDot
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYDotPtr() const
  {
    return yDot;
  }

  /** \fn void setYDot (const SimpleVector& newValue)
   *  \brief set the value of yDot to newValue
   *  \param SimpleVector newValue
   */
  inline void setYDot(const SimpleVector& newValue)
  {
    *yDot = newValue;
  }

  /** \fn void setYDotPtr(SimpleVector* newPtr)
   *  \brief set YDot to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setYDotPtr(SimpleVector *newPtr)
  {
    if (isYDotAllocatedIn) delete yDot;
    yDot = newPtr;
    isYDotAllocatedIn = false;
  }

  // -- lambda --

  /** \fn  const SimpleVector getLambda() const
   *  \brief get the value of lambda
   *  \return SimpleVector
   */
  inline const SimpleVector getLambda() const
  {
    return *lambda;
  }

  /** \fn SimpleVector* getLambdaPtr() const
   *  \brief get lambda
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLambdaPtr() const
  {
    return lambda;
  }

  /** \fn void setLambda (const SimpleVector& newValue)
   *  \brief set the value of lambda to newValue
   *  \param SimpleVector newValue
   */
  inline void setLambda(const SimpleVector& newValue)
  {
    *lambda = newValue;
  }

  /** \fn void setLambdaPtr(SimpleVector* newPtr)
   *  \brief set Lambda to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setLambdaPtr(SimpleVector *newPtr)
  {
    if (isLambdaAllocatedIn) delete lambda;
    lambda = newPtr;
    isLambdaAllocatedIn = false;
  }

  // -- yOld --

  /** \fn  const SimpleVector getYOld() const
   *  \brief get the value of yOld
   *  \return SimpleVector
   */
  inline const SimpleVector getYOld() const
  {
    return *yOld;
  }

  /** \fn SimpleVector* getYOldPtr() const
   *  \brief get yOld
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYOldPtr() const
  {
    return yOld;
  }

  /** \fn void setYOld (const SimpleVector& newValue)
   *  \brief set the value of yOld to newValue
   *  \param SimpleVector newValue
   */
  inline void setYOld(const SimpleVector& newValue)
  {
    *yOld = newValue;
  }

  /** \fn void setYOldPtr(SimpleVector* newPtr)
   *  \brief set YOld to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setYOldPtr(SimpleVector *newPtr)
  {
    if (isYOldAllocatedIn) delete yOld;
    yOld = newPtr;
    isYOldAllocatedIn = false;
  }

  // -- yDotOld --

  /** \fn  const SimpleVector getYDotOld() const
   *  \brief get the value of yDotOld
   *  \return SimpleVector
   */
  inline const SimpleVector getYDotOld() const
  {
    return *yDotOld;
  }

  /** \fn SimpleVector* getYDotOldPtr() const
   *  \brief get yDotOld
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getYDotOldPtr() const
  {
    return yDotOld;
  }

  /** \fn void setYDotOld (const SimpleVector& newValue)
   *  \brief set the value of yDotOld to newValue
   *  \param SimpleVector newValue
   */
  inline void setYDotOld(const SimpleVector& newValue)
  {
    *yDotOld = newValue;
  }

  /** \fn void setYDotOldPtr(SimpleVector* newPtr)
   *  \brief set YDotOld to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setYDotOldPtr(SimpleVector *newPtr)
  {
    if (isYDotOldAllocatedIn) delete yDotOld;
    yDotOld = newPtr;
    isYDotOldAllocatedIn = false;
  }

  // -- lambdaOld --

  /** \fn  const SimpleVector getLambdaOld() const
   *  \brief get the value of lambdaOld
   *  \return SimpleVector
   */
  inline const SimpleVector getLambdaOld() const
  {
    return *lambdaOld;
  }

  /** \fn SimpleVector* getLambdaOldPtr() const
   *  \brief get lambdaOld
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getLambdaOldPtr() const
  {
    return lambdaOld;
  }

  /** \fn void setLambdaOld (const SimpleVector& newValue)
   *  \brief set the value of lambdaOld to newValue
   *  \param SimpleVector newValue
   */
  inline void setLambdaOld(const SimpleVector& newValue)
  {
    *lambdaOld = newValue;
  }

  /** \fn void setLambdaOldPtr(SimpleVector* newPtr)
   *  \brief set LambdaOld to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setLambdaOldPtr(SimpleVector *newPtr)
  {
    if (isLambdaOldAllocatedIn) delete lambdaOld;
    lambdaOld = newPtr;
    isLambdaOldAllocatedIn = false;
  }

  /** \fn vector<int> getStatus(void)
   *  \brief get the status of this Interaction
   *  \return the value of status for this Interaction
   */
  inline const std::vector<int> getStatus() const
  {
    return status;
  }

  /** \fn void setStatus(vector<int>)
   *  \brief set the status of this Interaction
   *  \param the vector of integer to set the status
   */
  inline void setStatus(const std::vector<int>& vs)
  {
    status = vs;
  }

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

  /** \fn void check(const double& time, const double& pasH)
   * \brief compares the output of the relation with respect to the NonSmoothLaw and set the status
   * \param double : current time
   * \param double : current time step
   */
  void check(const double& time, const double& pasH);

  /** \fn void update(const double& time, const double& pasH)
   *  \brief set the status after the computation
   *  \param double : current time
   * \param double : current time step
   */
  void update(const double& time, const double& pasH);

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

private:
  /** \fn Interaction()
   *  \brief default constructor
   */
  Interaction();

  // --- MEMBERS ---
  /** name of the Interaction */
  std::string  id;
  /** number specific to each Interaction */
  int number;
  /** size of the the vector y */
  int nInteraction;
  /** relation between constrained variables and states variables */
  SimpleVector* y;
  /** relation between constrained variables and states variables */
  SimpleVector* yDot;
  /** result of the computeInput function */
  SimpleVector* lambda;
  /** relation between constrained variables and states variables */
  SimpleVector* yOld;
  /** relation between constrained variables and states variables */
  SimpleVector* yDotOld;
  /** result of the computeInput function */
  SimpleVector* lambdaOld;

  /** shows the status of the Interaction */
  std::vector<int> status;

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

  bool isYAllocatedIn;
  bool isYDotAllocatedIn;
  bool isLambdaAllocatedIn;
  bool isYOldAllocatedIn;
  bool isYDotOldAllocatedIn;
  bool isLambdaOldAllocatedIn;
  bool isRelationAllocatedIn;
  bool isNsLawAllocatedIn;
};

#endif // INTERACTION_H
