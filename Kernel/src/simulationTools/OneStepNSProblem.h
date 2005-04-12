#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "Strategy.h"
#include "Interaction.h"
#include "EqualityConstraint.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "OneStepNSProblemXML.h"

#include "SiconosConst.h"
#include "SiconosNumerics.h"

//using namespace std;

class Strategy;
class Interaction;
class EqualityConstraint;

class OneStepNSProblemXML;

extern string   DefaultSolver;
extern string   DefaultAlgoName;
extern string   DefaultAlgoNormType;
extern double   DefaultAlgoTolerance;
extern int    DefaultAlgoMaxIter;
extern double   DefaultAlgoSearchDirection;


/** \struct ConnectedInteraction
 *  \brief interaction connected to another interaction
 */
typedef struct
{
  /* state of the interaction connected :
   *  0->potential connection
   *  1->active connection
   */
  int status;

  /* the interaction connected to another interaction */
  Interaction *connected;

  /* position of the common DynamicalSystem in the dynamical system vector
   *  of the interaction in origin Interaction */
  int originInteractionDSRank;

  /* position of the common DynamicalSystem in the dynamical system vector
   *  of the interaction in connected Interaction */
  int connectedInteractionDSRank;
} Connection;

/** \class OneStepNSProblem
 *  \brief It's the part of the Strategy which solve the Interactions
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 */
class OneStepNSProblem
{
public:

  /** \fn OneStepNSProblem()
   *  \brief default constructor
   */
  OneStepNSProblem();

  /** \fn OneStepNSProblem(OneStepNSProblemXML*)
   *  \brief constructor with XML object of the OneStepNSProblem
   *  \param OneStepNSProblemXML* : the XML object corresponding
   */
  OneStepNSProblem(OneStepNSProblemXML*);

  virtual ~OneStepNSProblem();

  // getter/setter
  /** \fn int getN()
   *  \brief allow to get the value of n
   *  \return the value of n
   */
  inline int getN() const
  {
    return this->n;
  };

  /** \fn Strategy* getStrategy()
   *  \brief allow to get the Strategy
   *  \return the Strategy
   */
  inline Strategy* getStrategy() const
  {
    return this->strategy;
  }

  /** \fn vector< Interaction* > getInteractions()
   *  \brief allow to get the the vector of Interaction
   *  \return the vector interactionVector
   */
  inline vector< Interaction* > getInteractions() const
  {
    return this->interactionVector;
  };

  /** \fn Interaction* getInteraction(int)
   *  \brief allow to get a specific Interaction
   *  \param int the position of a specific Interaction in the vector of Interaction
   *  \return the specified Interaction
   */
  Interaction* getInteraction(const int);


  /** \fn void setN(int)
   *  \brief allow to set the value of n
   *  \param int : the value to set n
   */
  inline void setN(const int N)
  {
    this->n = N;
  };

  /** \fn void setInteractions(vector< Interaction* >)
   *  \brief allow to set the vector of Interactoin of the OneStepNSProblem
   *  \param vector<Interaction*> : the vector to set interactionVector
   */
  inline void setInteractions(const vector< Interaction* > interactions)
  {
    this->interactionVector = interactions;
  };

  /** \fn void addInteraction(Interaction*)
   *  \brief allow to add an Interaction to the OneStepNSProblem
   *  \param Interaction* : the Interaction to add to the vector of Interaction
   */
  void addInteraction(Interaction*);

  /** \fn inline OneStepNSProblemXML* getOneStepNSProblemXML()
   *  \brief allows to get the OneStepNSProblemXML of the OneStepNSProblem
   *  \return a pointer on the OneStepNSProblemXML of the OneStepNSProblem
   */
  inline OneStepNSProblemXML* getOneStepNSProblemXML() const
  {
    return this->onestepnspbxml;
  }

  /** \fn inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
   *  \brief allows to set the OneStepNSProblemXML of the OneStepNSProblem
   *  \param OneStepNSProblemXML* : the pointer to set OneStepNSProblemXML
   */
  inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
  {
    this->onestepnspbxml = osnspb;
  }

  /** \fn void setStrategy(Strategy*)
   *  \brief allow to set the Strategy of the OneStepNSProblem
   *  \return the Strategy*
   */
  inline void setStrategy(Strategy* str)
  {
    this->strategy = str;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the OneStepNSProblem
   *  \return string : the type of the OneStepNSProblem
   */
  inline string getType() const
  {
    return this->nspbType;
  }

  /////////////////////////////

  /** \fn void initialize(void)
  *  \brief initializes the OneStepNSProblem  and its interactions before doing the first step
  */
  virtual void initialize(void);

  /** \fn void nextStep(void)
  *  \brief prepares the Interaction for the next time step push y and lambda in Memory
  *  \exception to be defined
  *  \return void
  */
  virtual void nextStep(void);

  /** \fn void updateState(void)
   *  \brief predict all the relations before updating state of the problem
  */
  virtual void updateState(void);

  /** \fn void checkInteraction(void)
   *  \brief predict all the relations to see which ones have an effect
   */
  virtual void checkInteraction(void);

  /** \fn void formalize(void)
   *  \brief transform the discretised problem in a problem under numerical form
   *  param double : current time
   */
  virtual void formalize(double time);

  /** \fn void compute(void)
   *  \brief make the computation so solve the NS problem
   */
  virtual void compute(void);

  /** \fn void fillNSProblemWithNSProblemXML()
   *  \brief uses the OneStepNSProblemXML of the OneStepNSProblem to fill the fields of this OneStepNSProblem
   *  \exception RuntimeException
   */
  virtual void fillNSProblemWithNSProblemXML();

  /** \fn void fillSolvingMethod()
   *  \brief fills the fields of the solvingMethod object with data read in the XML DOM tree
   *  \exception RuntimeException
   */
  void fillSolvingMethod();

  /** \fn void saveNSProblemToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveNSProblemToXML();

  /** \fn void fillInteractionVector()
   *  \brief fill the interactionVector to store the Interactions concerned by this OneStepNSProblem
   *  This method is called when the platform is manually built and all the Interactions of the NonSmoothDynamicalSystem will be added to the interactionVector
   *  \exception RuntimeException
   */
  void fillInteractionVector();

  /** \fn void init(void)
   *  \brief initialize the value nInteraction if the input data don't do it. This value correspond to the size of the y vector
   */
  void init(void);

  /** \fn bool allInteractionConcerned()
  *   \brief defines if the vector of Interaction concerned by this OneStepNSProblem
  *          contains all the Interactions of the NonSmoothDynamicalSystem
  *   \return bool : true if the vector of Interaction of the OneStepNSProblem and the NonSmoothDynamicalSystem have the same size
  */
  bool allInteractionConcerned();

  /** \fn void setLemkeAlgorithm( double )
  *   \brief sets the parameters for the lemke solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the *tolerance* maxIter parameter
  */
  virtual void setLemkeAlgorithm(string, double = /*DefaultAlgoTolerance*/ DefaultAlgoMaxIter);

  /** \fn void setGsnlAlgorithm( double, string, int)
  *   \brief sets the parameters for the gsnl solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  */
  virtual void setGsnlAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
                                int = DefaultAlgoMaxIter);

  /** \fn void setGcpAlgorithm( double, string, int )
  *   \brief sets the parameters for the gcp solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  *   \param int : the iterMax parameter
  */
  virtual void setGcpAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
                               int = DefaultAlgoMaxIter);

  /** \fn void setLatinAlgorithm( double, string, int, double )
  *   \brief sets the parameters for the lcp solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  *   \param int : the iterMax parameter
  *   \param double : the search direction parameter
  */
  virtual void setLatinAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
                                 int = DefaultAlgoMaxIter, double = DefaultAlgoSearchDirection);


  /** \fn void updateConnectedInteractionMap()
  *   \brief mofifies the connectedInteractions map according to the interactions of the OneStepNSProblem
  */
  void updateConnectedInteractionMap();

  /** \fn void displayConnectedInteractionMap()
  *   \brief display the map of the connected interactions
  */
  void displayConnectedInteractionMap();


protected:
  /** type of the OneStepNSProblem */
  string nspbType;

  /** contains the size of the problem to solve */
  int n;

  Strategy *strategy;

  /** all the Interaction known by the OneStepNSProblem */
  vector < Interaction* > interactionVector;

  /** all the EqualityConstraint known by the OneStepNSProblem */
  vector < EqualityConstraint* > ecVector;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  OneStepNSProblemXML* onestepnspbxml;

  /** structure containing the structures of the numerous solving methods */
  methode solvingMethod;

  /** name of the solver to use */
  string solver;

  /** array of the connected interactions
   * in this map, we put all the active interactions (status = 1)
   * If an active interaction has no connection, the associated vector<Connection*> is empty */
  map< Interaction*, vector<Connection*> > connectedInteractionMap;
};

#endif // ONESTEPNSPROBLEM_H
