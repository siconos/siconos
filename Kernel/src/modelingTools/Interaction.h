#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>
#include <string>

#include "DynamicalSystem.h"
#include "NonSmoothLaw.h"
#include "Relation.h"

#include "InteractionXML.h"

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"

#include "SiconosConst.h"


//using namespace std;

class DynamicalSystem;
class Relation;
class NonSmoothLaw;

class InteractionXML;

/** \class Interaction
 *  \brief link between dynamic systems, driven by relations and non-smooth laws.
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 */
class Interaction
{
public:

  /** \fn Interaction()
   *  \brief default constructor : initialize pointers to NULL, integers to 0, etc.
   */
  Interaction();

  /** \fn Interaction(Interaction*)
   *  \brief copy constructor
   *  \param Interaction* : the object to copy
   */
  Interaction(Interaction *inter)
  {
    *this = *inter;
  };

  /** \fn Interaction(InteractionXML*)
   *  \brief constructor with XML object of the Interaction
   *  \param InteractionXML* : the XML object corresponding
   */
  Interaction(InteractionXML*);

  ~Interaction();

  // getter et setter

  /** \fn string getId(void)
   *  \brief allows to get the id of this Interaction
   *  \return the string, id of this Interaction
   */
  inline string getId(void) const
  {
    return this->id;
  };

  /** \fn int getNInteraction(void)
   *  \brief allows to get the number nInteraction of this Interaction
   *  \return the value of nInteraction
   */
  inline int getNInteraction(void) const
  {
    return this->nInteraction;
  };

  /** \fn SimpleVector getY(void)
   *  \brief allows to get the SiconosVector y of this Interaction
   *  \return SimpleVector, y of this Interaction
   */
  inline SimpleVector getY(void) const
  {
    return this->y;
  };

  /** \fn SiconosVector* getYPtr(void)
   *  \brief allows to get the SiconosVector* y of this Interaction
   *  \return SiconosVector*, y of this Interaction
   */
  SimpleVector* getYPtr(void);

  /** \fn SimpleVector getLambda(void)
   *  \brief allows to get the SiconosVector lambda of this Interaction
   *  \return SimpleVector, lambda of this Interaction
   */
  inline SimpleVector getLambda(void) const
  {
    return this->lambda;
  };

  /** \fn SiconosVector* getLambdaPtr(void)
   *  \brief allows to get the SiconosVector* lambda of this Interaction
   *  \return SiconosVector*, lambda of this Interaction
   */
  SimpleVector* getLambdaPtr(void);

  /** \fn SimpleVector getYOld(void)
   *  \brief allows to get the SiconosVector yOld of this Interaction
   *  \return SimpleVector, yOld of this Interaction
   */
  inline SimpleVector getYOld(void) const
  {
    return this->yOld;
  };

  /** \fn SimpleVector* getYOldPtr(void)
   *  \brief allows to get the SimpleVector* yOld of this Interaction
   *  \return SimpleVector*, yOld of this Interaction
   */
  SimpleVector* getYOldPtr(void);

  /** \fn inline void setYDot(SimpleVector &yDot)
   *  \brief set yDot of this Interaction
   *  \param SimpleVector &yDot : new value of yDot
   */
  inline void setYDot(const SimpleVector &yDot)
  {
    this->yDot = yDot;
  };

  /** \fn inline SimpleVector getYDot(void)
   *  \brief get yDot of this Interaction
   *  \return SimpleVector, yDot of this Interaction
   */
  inline  SimpleVector getYDot(void) const
  {
    return this->yDot;
  };

  /** \fn inline SimpleVector* getYDotPtr(void)
   *  \brief get yDot of this Interaction
   *  \return SimpleVector*, pointer on yDot of this Interaction
   */
  inline SimpleVector* getYDotPtr(void)
  {
    return &(this->yDot);
  };

  /** \fn inline void setYDotOld(SimpleVector &yDot)
   *  \brief set previous value of yDot of this Interaction
   *  \param SimpleVector &yDotOld : previous value of yDot
   */
  inline void setYDotOld(const SimpleVector& yDotOld)
  {
    this->yDotOld = yDotOld;
  };

  /** \fn inline SimpleVector getYDotOld(void)
   *  \brief get previous value of yDot of this Interaction
   *  \return SimpleVector, previous value of yDot of this Interaction
   */
  inline  SimpleVector getYDotOld(void) const
  {
    return this->yDotOld;
  };

  /** \fn inline SimpleVector* getYDotOldPtr(void)
   *  \brief get previous value of yDot of this Interaction
   *  \return SimpleVector*, pointer on previous value of yDot of this Interaction
   */
  inline SimpleVector* getYDotOldPtr(void)
  {
    return &(this->yDotOld);
  };

  /** \fn SimpleVector getLambdaOld(void)
   *  \brief get the SimpleVector lambdaOld of this Interaction
   *  \return SimpleVector, lambdaOld of this Interaction
   */
  inline SimpleVector getLambdaOld(void) const
  {
    return this->lambdaOld;
  };

  /** \fn SimpleVector* getLambdaOldPtr(void)
   *  \brief get the SimpleVector* lambdaOld of this Interaction
   *  \return SimpleVector*, lambdaOld of this Interaction
   */
  SimpleVector* getLambdaOldPtr(void);

  /** \fn vector<int> getStatus(void)
   *  \brief allows to get the status of this Interaction
   *  \return the value of status for this Interaction
   */
  inline vector<int> getStatus(void) const
  {
    return this->status;
  };

  /** \fn vector<DynamicalSystem*> getDynamicalSystems(void)
   *  \brief allows to get the 2 DynamicalSystem of this Interaction
   *  \return vector<DynamicalSystem*> : a vector of DS, with 2 DynamicalSystem
   */
  inline vector<DynamicalSystem*> getDynamicalSystems(void)
  {
    return this->vectorDS;
  };

  /** \fn DynamicalSystem* getDynamicalSystem(int)
   *  \brief allows to get a specific DynamicalSystem
   *  \param the identification number of the wanted DynamicalSystem
   *  \return null if there's no DS corresponding else the right DynamicalSystem
   */
  DynamicalSystem* getDynamicalSystem(int);

  /** \fn Relation* getRelation(void)
   *  \brief allows to get the Relation of this Interaction
   *  \return a pointer on this Relation
   */
  inline Relation* getRelation(void) const
  {
    return this->relation;
  };

  /** \fn NonSmoothLaw* getNonSmoothLaw(void)
   *  \brief allows to get the NonSmoothLaw of this Interaction
   *  \return a pointer on this NonSmoothLaw
   */
  inline NonSmoothLaw* getNonSmoothLaw(void) const
  {
    return this->nslaw;
  };

  /** \fn int getNumber()
   *  \brief allows to get the value of number
   *  \return the value of number
   */
  inline int getNumber() const
  {
    return this->number;
  };

  /** \fn void setNumber(int)
   *  \brief allows to modify the value of number
   *  \param int : the value to set for number
   */
  inline void setNumber(const int nb)
  {
    this->number = nb;
  };

  /** \fn void setId(int)
   *  \brief allows to set the id of this Interaction
   *  \param the integer to set the id
   */
  inline void setId(const int id)
  {
    this->id = id;
  };

  /** \fn void setNInteraction(int)
   *  \brief allows to set the value of nInteraction
   *  \param the integer to set the nInteraction
   */
  inline void setNInteraction(const int nInter)
  {
    this->nInteraction = nInter;
  };

  /** \fn void setLambda(SimpleVector)
   *  \brief allows to set the SimpleVector lambda
   *  \param SimpleVector : new value of lambda
   */
  inline void setLambda(const SimpleVector& lambda)
  {
    this->lambda = lambda;
  };

  /** \fn void setStatus(vector<int>)
   *  \brief allows to set the status of this Interaction
   *  \param the vector of integer to set the status
   */
  inline void setStatus(vector<int> vs)
  {
    this->status = vs;
  };

  /** \fn void setDynamicalSystems(DynamicalSystem*, DynamicalSystem*)
   *  \brief allows to set the 2 DS of this Interaction
   *  \param the first DS* to set in the Interaction
   *  \param the second DS* to set in the Interaction
   */
  void setDynamicalSystems(DynamicalSystem*, DynamicalSystem*);

  /** \fn void setRelation(Relation*)
   *  \brief allows to set the Relation of this Interaction
   *  \param the Relation* to set
   */
  inline void setRelation(Relation*)
  {
    this->relation = relation;
  };

  /** \fn void setNonSmoothLaw(NonSmoothLaw*)
   *  \brief allows to set the NonSmoothLaw of this Interaction
   *  \param the NonSmoothLaw* to set
   */
  inline void setNonSmoothLaw(NonSmoothLaw*)
  {
    this->nslaw = nslaw;
  };

  /** \fn inline InteractionXML* getInteractionXML()
   *  \brief allows to get the InteractionXML* of the Interaction
   *  \return InteractionXML* : the pointer on the InteractionXML
   */
  inline InteractionXML* getInteractionXML()
  {
    return this->interactionxml;
  }

  /** \fn inline void setInteractionXML(InteractionXML* interxml)
   *  \brief allows to set the InteractionXML* of the Interaction
   *  \param InteractionXML* :  the pointer to set
   */
  inline void setInteractionXML(InteractionXML* interxml)
  {
    this->interactionxml = interxml;
  }


  ////////////////////////////////

  /** \fn   void initialize(void);
   * \brief   initialization of the interaction
   * \warning : this function does nothing apart displaying data of the Interaction
   */
  void initialize();

  /** \fn   void swapInMemory(void);
   * \brief   put values of y in yOld, the same for lambda
   */
  void swapInMemory(void);

  /** \fn void check(double time)
   * \brief compares the output of the relation with respect to the NonSmoothLaw and set the status
   * \param double : current time
   */
  void check(double time);

  /** \fn void update(double time)
   *  \brief set the status after the computation
   *  \param double : current time
   */
  void update(double time);

  /** \fn void saveInteractionToXML()
   *  \brief copy the data of the Interaction to the XML tree
   *  \exception RuntimeException
   */
  void saveInteractionToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;


  /** \fn void createInteraction(InteractionXML * interactionXML, int number, int nInter,
              vector<int>* status, vector<DynamicalSystem*>* dsConcerned)
   *  \brief allows to create the Interaction with an xml file, or the needed data
   *  \param InteractionXML * : the XML object for this Interaction
   *  \param int : the number of this Interaction
   *  \param int : the value of nInter of this Interaction
   *  \param vector<int>* : the status of this Interaction
   *  \param vectir<DynamicalSystem*>* : the Dynamical Systems concerned by this interaction
   *  \exception RuntimeException
   */
  void createInteraction(InteractionXML * interactionXML, int number = -1, int nInter = -1,
                         vector<int>* status = NULL, vector<DynamicalSystem*>* dsConcerned = NULL); //, NonSmoothDynamicalSystem * nsds = NULL);

  /** \fn Relation* createLagrangianLinearR(SiconosMatrix* H, SiconosVector* b)
   *  \brief allows to create a LagrangianLinearR relation for this Interaction
   *  \param SiconosMatrix* : the H matrix
   *  \param SiconosVector* : the b vector
   *  \return Relation* : the relation created
   */
  Relation* createLagrangianLinearR(SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** \fn Relation* createLagrangianNonLinearR()
   *  \brief allows to create a LagrangianNonLinearR relation for this Interaction
   *  \param string : the name of the computeInput plugin
   *  \param string : the name of the computeOutput plugin
   *  \return Relation* : the relation created
   */
  Relation* createLagrangianNonLinearR(string computeInput, string computeOutput);

  /** \fn Relation* createLinearTIR(SiconosMatrix* C, SiconosMatrix* D,
                    SiconosMatrix* E, SiconosVector* a)
   *  \brief allows to create a LinearTIR relation for this Interaction
   *  \param SiconosMatrix* : the C matrix
   *  \param SiconosMatrix* : the D matrix
   *  \param SiconosMatrix* : the E matrix
   *  \param SiconosVector* : the a vector
   *  \return Relation* : the relation created
   */
  Relation* createLinearTIR(SiconosMatrix* C = NULL, SiconosMatrix* D = NULL,
                            SiconosMatrix* E = NULL, SiconosVector* a = NULL);


  /** \fn NonSmoothLaw* createComplementarityConditionNSL()
   *  \brief allows to create a ComplementarityConditionNSL non-smooth law for this Interaction
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createComplementarityConditionNSL();

  /** \fn NonSmoothLaw* createRelayNSL(double c, double d)
   *  \brief allows to create a RelayNSL non-smooth law for this Interaction
   *  \param double : the c value
   *  \param double : the d value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createRelayNSL(double c, double d);

  /** \fn NonSmoothLaw* createNewtonImpactLawNSL(double e)
   *  \brief allows to create a NewtonImpactLawNSL non-smooth law for this Interaction
   *  \param double : the e value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createNewtonImpactLawNSL(double e);

  /** \fn NonSmoothLaw* createNewtonImpactLawNSL(double e)
   *  \brief allows to create a NewtonImpactLawNSL non-smooth law for this Interaction
   *  \param double : the en value
   *  \param double : the et value
   *  \param double : the mu value
   *  \return NonSmoothLaw* : the non-smooth law created
   */
  NonSmoothLaw* createNewtonImpactFrictionNSL(double en, double et, double mu);


protected:
  /** \fn void fillInteractionWithInteractionXML()
   *  \brief uses the InteractionXML of the Interaction to fill the fields of this Interaction
   *  \exception RuntimeException
   */
  void fillInteractionWithInteractionXML();

  /** \fn void linkInteractionWithInteractionXML()
   *  \brief makes the links between the RelationXMLs, NonSmoothLawXMLs of the InteractionXML of the Interaction and the Relations, NonSmoothLaws
   */
  void linkInteractionWithInteractionXML();


private:
  /** name of the Interaction */
  string id;
  /** unique number specific to each Interaction */
  int number;
  /** size of the the vector y */
  int nInteraction;
  /** relation between constrained variables and states variables */
  SimpleVector y;
  /** result of the computeInput function */
  SimpleVector lambda;
  /** relation between constrained variables and states variables */
  SimpleVector yOld;
  /** relation between constrained variables and states variables */
  SimpleVector yDot;
  /** relation between constrained variables and states variables */
  SimpleVector yDotOld;
  /** result of the computeInput function */
  SimpleVector lambdaOld;


  /** shows the status of the Interaction */
  vector<int> status;

  /** the Dynamical Systems concerned by this relation
   * dynamical systems are given by pair :
   * the vector (ds1, ds2, ds3, ds4, ds5, ds6)
   * means ds1 interacts with ds2, ds3 interacts with ds4, ds5 interacts with ds6*/
  vector<DynamicalSystem*> vectorDS;

  /** the Non-smooth Law of the interaction*/
  NonSmoothLaw *nslaw;

  /** the type of Relation of the interaction */
  Relation *relation;

  /** the XML object linked to the Interaction to read XML data */
  InteractionXML *interactionxml;

};

#endif // INTERACTION_H
