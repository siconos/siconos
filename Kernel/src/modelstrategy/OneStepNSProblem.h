//$Id: OneStepNSProblem.h,v 1.41 2005/03/23 15:03:55 jbarbier Exp $
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

using namespace std;

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
 *  \author Jean-Michel Barbier
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 * $Date: 2005/03/23 15:03:55 $
 * $Revision: 1.41 $
 * $Author: jbarbier $
 * $Source: /CVS/Siconos/SICONOS/src/modelstrategy/OneStepNSProblem.h,v $
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
  void setLemkeAlgorithm(string, double = /*DefaultAlgoTolerance*/ DefaultAlgoMaxIter);

  /** \fn void setGsnlAlgorithm( double, string, int)
  *   \brief sets the parameters for the gsnl solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  */
  void setGsnlAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
                        int = DefaultAlgoMaxIter);

  /** \fn void setGcpAlgorithm( double, string, int )
  *   \brief sets the parameters for the gcp solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  *   \param int : the iterMax parameter
  */
  void setGcpAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
                       int = DefaultAlgoMaxIter);

  /** \fn void setLatinAlgorithm( double, string, int, double )
  *   \brief sets the parameters for the lcp solfing algorithm in the solvingMethod structure
  *   \param string : the kind of solving method
  *   \param double : the tolerance parameter
  *   \param string : the norm type paramater
  *   \param int : the iterMax parameter
  *   \param double : the search direction parameter
  */
  void setLatinAlgorithm(string, double = DefaultAlgoTolerance, string = DefaultAlgoNormType,
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
//$Log: OneStepNSProblem.h,v $
//Revision 1.41  2005/03/23 15:03:55  jbarbier
//- adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//
//Revision 1.40  2005/03/21 16:48:03  jbarbier
//- EqualityConstraint : computeInput and computeOutput functions added (plugin funcitons)
//
//- link OneStepNSProblem - EqualityConstraint established
//
//- modification of OneStepNSProblem save according to change to normType[64] in SiconosNumerics.h
//
//Revision 1.39  2005/03/21 10:17:11  jbarbier
//- normType available after the modification in Numerics
//
//- searchDirection save fixed when default value was used
//
//Revision 1.38  2005/03/02 16:06:34  jbarbier
//- DoubleContact sample runnig successfully!
//
//- computeM and computeQ of LCP fixed
//
//Revision 1.37  2005/02/14 09:52:21  charlety
//_ getters / setters put inline
//
//Revision 1.36  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.35  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.34  2005/02/01 10:43:48  jbarbier
//- BouncingBall sample taking account of the solver defined in the XML input file
//
//Revision 1.33  2005/01/25 09:27:18  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.32  2005/01/24 14:33:02  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.31  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.30  2005/01/14 11:51:25  jbarbier
//- attribute "all" of the OneStepNSProblem terminated
//in OneStepIntegrator and Interaction the attribute is available
//
//Revision 1.29  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.28  2004/09/22 14:11:14  charlety
//
//  _ revision of Doxygen comments in modelstrategy
//
//Revision 1.27  2004/09/21 11:49:10  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.26  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.25  2004/09/14 13:49:54  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.24  2004/09/10 11:26:17  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.23  2004/08/11 14:43:45  jbarbier
//- beginning of the mechanism of creation without XML input file of the objects of the platform with the
//creatObjtect methods
//
//- function saveWToXML for Moreau integrator, and same specific functions to save
//M,q and Q,p for LCP and QP
//
//- function to check coherency of the Model
//
//Revision 1.22  2004/07/29 14:25:40  jbarbier
//- $Log: OneStepNSProblem.h,v $
//- Revision 1.41  2005/03/23 15:03:55  jbarbier
//- - adaptation to the LMGC90 tags in non smooth dynamical system and strategy
//-
//- Revision 1.40  2005/03/21 16:48:03  jbarbier
//- - EqualityConstraint : computeInput and computeOutput functions added (plugin funcitons)
//-
//- - link OneStepNSProblem - EqualityConstraint established
//-
//- - modification of OneStepNSProblem save according to change to normType[64] in SiconosNumerics.h
//-
//- Revision 1.39  2005/03/21 10:17:11  jbarbier
//- - normType available after the modification in Numerics
//-
//- - searchDirection save fixed when default value was used
//-
//- Revision 1.38  2005/03/02 16:06:34  jbarbier
//- - DoubleContact sample runnig successfully!
//-
//- - computeM and computeQ of LCP fixed
//-
//- Revision 1.37  2005/02/14 09:52:21  charlety
//- _ getters / setters put inline
//-
//- Revision 1.36  2005/02/10 10:35:19  jbarbier
//- - new file regrouping all the const values of the model, modelingTools and numericalStrategy
//-
//- - new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//-
//- - new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//-
//- Revision 1.35  2005/02/04 14:52:44  jbarbier
//- - Rolling balls in progress (contact is detected)
//-
//- - time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//-
//- Revision 1.34  2005/02/01 10:43:48  jbarbier
//- - BouncingBall sample taking account of the solver defined in the XML input file
//-
//- Revision 1.33  2005/01/25 09:27:18  jbarbier
//- - save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//-
//- Revision 1.32  2005/01/24 14:33:02  jbarbier
//- - OneStepNSProblem > Solver tag is available and managed in the XML part
//-
//- - tests added on OneStepNSProblem > Solver tag
//-
//- Revision 1.31  2005/01/20 14:44:49  jbarbier
//- - NSDS class renamed NonSmoothDynamicalSystem
//-
//- - code reduce, some comments remove
//-
//- Revision 1.30  2005/01/14 11:51:25  jbarbier
//- - attribute "all" of the OneStepNSProblem terminated
//- in OneStepIntegrator and Interaction the attribute is available
//-
//- Revision 1.29  2004/12/06 10:10:34  jbarbier
//- - integration of Numerics and use of Numerics on the bouncing ball sample
//-
//- - Numerics is now needed to run the bouncing ball sample!
//-
//- Revision 1.28  2004/09/22 14:11:14  charlety
//-
//-   _ revision of Doxygen comments in modelstrategy
//-
//- Revision 1.27  2004/09/21 11:49:10  jbarbier
//- - correction in the XML save for a manual construction of the platform :
//-     DS_Concerned of the Interaction
//-     DS_Concerned of the Integrator
//-
//- - test updated for these changes
//-
//- Revision 1.26  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.25  2004/09/14 13:49:54  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.24  2004/09/10 11:26:17  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.23  2004/08/11 14:43:45  jbarbier
//- - beginning of the mechanism of creation without XML input file of the objects of the platform with the
//- creatObjtect methods
//-
//- - function saveWToXML for Moreau integrator, and same specific functions to save
//- M,q and Q,p for LCP and QP
//-
//- - function to check coherency of the Model
//- and $Id: OneStepNSProblem.h,v 1.41 2005/03/23 15:03:55 jbarbier Exp $ added
//
