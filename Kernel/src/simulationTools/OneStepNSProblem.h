#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "Strategy.h"
#include "Interaction.h"
#include "EqualityConstraint.h"
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "OneStepNSProblemXML.h"
#include "Topology.h"
#include "SiconosConst.h"
#include "SiconosNumerics.h"
#include "check.h"
#include <iostream>
#include <vector>

class Strategy;
class Interaction;
class EqualityConstraint;

class OneStepNSProblemXML;

extern std::string  DefaultSolver;
extern std::string  DefaultAlgoName;
extern std::string  DefaultAlgoNormType;
extern double   DefaultAlgoTolerance;
extern int  DefaultAlgoMaxIter;
extern double   DefaultAlgoSearchDirection;

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

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn OneStepNSProblem()
   *  \brief default constructor
   */
  OneStepNSProblem();

  /** \fn OneStepNSProblem(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  OneStepNSProblem(OneStepNSProblemXML*, Strategy * = NULL);

  /** \fn OneStepNSProblem(Strategy*, const string& solverName, const string& newSolvingMethod, const int& maxIter,
   *                       const double & Tolerance=0, const string & NormType="none",
   *                       const double & SearchDirection=0)
   *  \brief constructor from data
   *  \param Strategy *: the strategy that owns this problem
   *  \param string: solver name (optional)
   *  \param string: name of the solving method (optional but required if a solver is given)
   *  \param int : MaxIter (optional) required if a solver is given
   *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin
   *  \param string : NormType (optional) -> for NLGS, Gcp, Latin
   *  \param double : SearchDirection (optional) -> for Latin
   */
  OneStepNSProblem(Strategy * , const std::string& = "none", const std::string& = "none",
                   const int& = 0, const double& = 0, const std::string & = "none",
                   const double & = 0);

  virtual ~OneStepNSProblem();

  // --- GETTERS/SETTERS ---

  /** \fn inline const string getType() const
   *  \brief to get the type of the OneStepNSProblem
   *  \return string
   */
  inline std::string getType() const
  {
    return nspbType;
  }

  /** \fn inline void setType(const string&)
   *  \brief set the type of the OneStepNSProblem
   *  \param: string
   */
  inline void setType(const std::string & newVal)
  {
    nspbType = newVal;
  }

  /** \fn const int getN() const
   *  \brief get the value of n
   *  \return an int
   */
  inline const unsigned int getN() const
  {
    return n;
  }

  /** \fn void setN(const int&)
   *  \brief set the value of n
   *  \param an int
   */
  inline void setN(const unsigned int& newVal)
  {
    n = newVal;
  }

  /** \fn vector< Interaction* > getInteractions()
   *  \brief get the the vector of Interaction
   *  \return a vector stl
   */
  inline const std::vector< Interaction* > getInteractions() const
  {
    return interactionVector;
  }

  /** \fn Interaction* getInteractionPtr(const int&)
   *  \brief get a specific Interaction
   *  \param int the position of a specific Interaction in the vector of Interaction
   *  \return a pointer on Interaction
   */
  Interaction* getInteractionPtr(const unsigned int&);

  /** \fn void setInteractions(vector< Interaction* >)
    *  \brief set the vector of Interaction
    *  \param vector<Interaction*> : a vector stl
    */
  inline void setInteractions(const std::vector< Interaction* >& newVec)
  {
    interactionVector = newVec;
  }

  /** \fn vector< EqualityConstraint* > getEqualityConstraints()
   *  \brief get the the vector of EqualityConstraint
   *  \return a vector stl
   */
  inline const std::vector< EqualityConstraint* > getEqualityConstraints() const
  {
    return ecVector;
  }

  /** \fn void setEqualityConstraints(vector< EqualityConstraint* >)
    *  \brief set the vector of EqualityConstraint
    *  \param vector<EqualityConstraint*> : a vector stl
    */
  inline void setEqualityConstraints(const std::vector< EqualityConstraint* >& newVec)
  {
    ecVector = newVec;
  }

  /** \fn inline const string getSolver() const
   *  \brief to get the solver of the OneStepNSProblem
   *  \return string
   */
  inline std::string getSolver() const
  {
    return solver;
  }

  /** \fn inline void setSolver(const string&)
   *  \brief set the solver of the OneStepNSProblem
   *  \param: string
   */
  inline void setSolver(const std::string & newVal)
  {
    solver = newVal;
  }

  /** \fn Strategy* getStrategyPtr()
   *  \brief get the Strategy
   *  \return a pointer on Strategy
   */
  inline Strategy* getStrategyPtr() const
  {
    return strategy;
  }

  /** \fn void setStrategyPtr(Strategy*)
   *  \brief set the Strategy of the OneStepNSProblem
   *  \param: a pointer on Strategy
   */
  inline void setStrategy(Strategy* str)
  {
    strategy = str;
  }

  /** \fn inline OneStepNSProblemXML* getOneStepNSProblemXML()
   *  \brief get the OneStepNSProblemXML
   *  \return a pointer on OneStepNSProblemXML
   */
  inline OneStepNSProblemXML* getOneStepNSProblemXML() const
  {
    return onestepnspbxml;
  }

  /** \fn inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
   *  \brief set the OneStepNSProblemXML
   *  \param a pointer on OneStepNSProblemXML
   */
  inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
  {
    onestepnspbxml = osnspb;
  }

  // --- OTHER FUNCTIONS ---

  /** \fn void addInteraction(Interaction*)
   *  \brief add an Interaction to the OneStepNSProblem
   *  \param Interaction* : the Interaction to add to the vector of Interaction
   */
  void addInteraction(Interaction*);

  /** \fn void initialize()
   *  \brief initialize the problem(compute topology ...)
   */
  virtual void initialize();

  /** \fn void computeEffectiveOutput();
   *  \brief compute variables indexMax, effectiveOutput and effectiveSizeOutput of the topology of the nsds
   */
  void computeEffectiveOutput();

  /** \fn void nextStep(void)
   *  \brief prepares the problem for the next time step
   *  \exception to be defined
   */
  virtual void nextStep();

  /** \fn void updateInput()
   *  \brief compute r thanks to lambda
   */
  virtual void updateInput();

  /** \fn void updateOutput(void)
   *  \brief compute output for all the interactions
   */
  virtual void updateOutput();

  /** \fn void compute(const double&)
   *  \brief make the computation so solve the NS problem
   *  param double : current time
   */
  virtual void compute(const double&);

  /** \fn void fillSolvingMethod(const string& newSolvingMethod, const int& maxIter,
   *                             const double & Tolerance=0, const string & NormType="none",
   *                             const double & SearchDirection=0)
   *  \brief to fill the fields of the solvingMethod object
   *  \param string: name of the solving method
   *  \param int : MaxIter
   *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin
   *  \param string : NormType (optional) -> for NLGS, Gcp, Latin
   *  \param double : SearchDirection (optional) -> for Latin
   *  \exception RuntimeException
   */
  void fillSolvingMethod(const std::string&, const int& ,
                         const double& = 0, const std::string & = "none",
                         const double & = 0);

  /** \fn void saveNSProblemToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveNSProblemToXML();

  /** \fn bool allInteractionConcerned()
   *   \brief defines if the vector of Interaction concerned by this OneStepNSProblem
   *          contains all the Interactions of the NonSmoothDynamicalSystem
   *   \return bool : true if the vector of Interaction of the OneStepNSProblem and the NonSmoothDynamicalSystem have the same size
   */
  bool allInteractionConcerned();

  /** \fn void setLemkeAlgorithm(const std::string&, const double& = DefaultAlgoMaxIter)
   *   \brief sets the parameters for the lemke solving algorithm in the solvingMethod structure
   *   \param string : the solving method
   *   \param unsigned int : the max number of iterations
   */
  virtual void setLemkeAlgorithm(const std::string&, const unsigned int& = DefaultAlgoMaxIter);

  /** \fn void setLexicoLemkeAlgorithm(const std::string&, const double& = DefaultAlgoMaxIter)
   *   \brief sets the parameters for the lexicoLemke solving algorithm in the solvingMethod structure
   *   \param string : the solving method
   *   \param  unsigned int: the max number of iterations
   */
  virtual void setLexicoLemkeAlgorithm(const std::string&, const unsigned int& = DefaultAlgoMaxIter);

  /** \fn void setNLGSAlgorithm(const std::string&, const double& = DefaultAlgoTolerance, const std::string& = DefaultAlgoNormType,
   *                            const int& = DefaultAlgoMaxIter)
   *   \brief sets the parameters for the NLGS solving algorithm in the solvingMethod structure
   *   \param string : solving method
   *   \param double : tolerance
   *   \param int    : max number of iterations
   *   \param string : norm type (unused)
   */
  virtual void setNLGSAlgorithm(const std::string&, const double& = DefaultAlgoTolerance, const int& = DefaultAlgoMaxIter,
                                const std::string& = DefaultAlgoNormType);

  /** \fn void setQPAlgorithm(const std::string&, const double& = DefaultAlgoTolerance)
   *   \brief sets the parameters for the QP solving algorithm in the solvingMethod structure
   *   \param string : solving method
   *   \param double : tolerance
   */
  virtual void setQPAlgorithm(const std::string&, const double& = DefaultAlgoTolerance);

  /** \fn void setNSQPAlgorithm(const std::string&, const double& = DefaultAlgoTolerance)
   *   \brief sets the parameters for the NSQP solving algorithm in the solvingMethod structure
   *   \param string : solving method
   *   \param double : tolerance
   */
  virtual void setNSQPAlgorithm(const std::string&, const double& = DefaultAlgoTolerance);

  /** \fn void setCPGAlgorithm( const std::string&, const double& = DefaultAlgoTolerance, const std::string& = DefaultAlgoNormType,
   *                            const int& = DefaultAlgoMaxIter )
   *   \brief sets the parameters for the gcp solving algorithm in the solvingMethod structure
   *   \param string : solving method
   *   \param double : tolerance
   *   \param int    : max number of iterations
   *   \param string : norm type
   */
  virtual void setCPGAlgorithm(const std::string&, const double& = DefaultAlgoTolerance, const int& = DefaultAlgoMaxIter ,
                               const std::string& = DefaultAlgoNormType);

  /** \fn void setLatinAlgorithm(const std::string&, const double& = DefaultAlgoTolerance, const int& = DefaultAlgoMaxIter,
   *                             const std::string& = DefaultAlgoNormTypeconst, double& = DefaultAlgoSearchDirection)
   *   \brief sets the parameters for the lcp solfing algorithm in the solvingMethod structure
   *   \param string : solving method
   *   \param double : tolerance
   *   \param int    : iterMax
   *   \param string : norm type (unused)
   *   \param double : the search direction parameter
   */
  virtual void setLatinAlgorithm(const std::string&, const double& = DefaultAlgoTolerance, const int& = DefaultAlgoMaxIter,
                                 const std::string& = DefaultAlgoNormType, const double& = DefaultAlgoSearchDirection);

  bool isOneStepNsProblemComplete();

protected:
  /** type of the OneStepNSProblem */
  std::string nspbType;

  /** size of the problem to solve */
  unsigned int n;

  /** all the Interaction known by the OneStepNSProblem */
  std::vector< Interaction* > interactionVector;

  /** all the EqualityConstraint known by the OneStepNSProblem */
  std::vector< EqualityConstraint* > ecVector;

  /** map that links each interaction with the corresponding diagonal block*/
  std::map< Interaction* , SiconosMatrix*>  diagonalBlocksMap  ;

  /** map that links each interaction with the corresponding extra-diagonal blocks
      map < InteractionA * , map < InteractionB* , blockMatrixAB > >
      Interaction A and B are coupled through blockMatrixAB.
  */
  std::map< Interaction* , std::map<Interaction *, SiconosMatrix*> >  extraDiagonalBlocksMap  ;

  /** structure containing the structures of the numerous solving methods */
  method solvingMethod;

  /** name of the solver to use */
  std::string solver;

  /** link to the strategy that owns the NSPb */
  Strategy *strategy;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  OneStepNSProblemXML* onestepnspbxml;
};

#endif // ONESTEPNSPROBLEM_H
