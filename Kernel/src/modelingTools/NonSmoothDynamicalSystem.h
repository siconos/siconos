#ifndef NSDS_H
#define NSDS_H

#include <iostream>
#include <vector>
#include <string>
#include <typeinfo>
#include <libxml/tree.h>

#include "Model.h"
#include "DynamicalSystem.h"
#include "BoundaryCondition.h"
#include "Interaction.h"
#include "EqualityConstraint.h"
#include "NSDSXML.h"

#include "SiconosConst.h"

//using namespace std;


enum dynamicalsystem {LAGRANGIANNLDS, LAGRANGIANTIDS, LINEARTIDS};

class Model;
class Interaction;
class DynamicalSystem;
class BoundaryCondition;
class EqualityConstraint;
class NSDSXML;

/** \class NonSmoothDynamicalSystem
 *  \brief This class contains the formalization of the Non Smooth Dynamical System (NonSmoothDynamicalSystem)
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) Apr 23, 2004
 *
 *
 *
 *   \warning
 */
class NonSmoothDynamicalSystem
{
public:

  /** \fn NonSmoothDynamicalSystem()
   *  \brief default constructor
   *  \param NSDSXML* : the XML object corresponding to the NonSmoothDynamicalSystem
   */
  NonSmoothDynamicalSystem();

  /** \fn NonSmoothDynamicalSystem(NSDSXML*)
   *  \brief constructor with XML object of the NonSmoothDynamicalSystem
   *  \param  const NonSmoothDynamicalSystem* : the XML object corresponding to the NonSmoothDynamicalSystem
   */
  NonSmoothDynamicalSystem(NSDSXML*);

  /** \fn NonSmoothDynamicalSystem(bool)
   *  \brief constructor with indication concerning boundary conditions
   *  \param const bool : if true, the systems has boundary conditions
   */
  NonSmoothDynamicalSystem(bool);

  /** \fn NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem*)
   *  \brief constructor by copy
   *  \param  const NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem object which must be copied
   */
  NonSmoothDynamicalSystem(NonSmoothDynamicalSystem*);

  /** \fn NonSmoothDynamicalSystem(string)
   *  \brief constructs the object with is type (IVP or BVP)
   *  \param a string value which determines if the problem is BVP or IVP
   */
  NonSmoothDynamicalSystem(string);

  ~NonSmoothDynamicalSystem();

  /** \fn bool isBVP(void)
   *  \brief determines if the NonSmoothDynamicalSystem is a Boundary Value Problem (BVP)
   *  \return the value of BVP
   */
  inline bool isBVP(void) const
  {
    return this->BVP;
  };

  /** \fn bool isIVP(void)
   *  \brief determines if the NonSmoothDynamicalSystem is an Initial Value Problem (IVP)
   *  \return the value of IVP
   */
  inline bool isIVP(void) const
  {
    return !this->BVP;
  };

  /** \fn inline int getDSVectorSize()
   *  \brief allows to get the size of the vector of Dynamical Systems
   *  \return int : the size of DSVector
   */
  inline int getDSVectorSize() const
  {
    return this->DSVector.size();
  }

  /** \fn vector<DynamicalSystem*> getDynamicalSystems(void)
   *  \brief allows to get all the DynamicalSystem of the NonSmoothDynamicalSystem problem
   *  \return the vector of DS
   */
  inline vector<DynamicalSystem*> getDynamicalSystems(void) const
  {
    return this->DSVector;
  };

  /** \fn DynamicalSystem* getDynamicalSystem(int)
   *  \brief allows to get one specific DynamicalSystem, with its place in the DSVector
   *  \param the number in the vector of DynamicalSystem to get
   *  \return the DynamicalSystem at the right place in the DSVector
   */
  DynamicalSystem* getDynamicalSystem(int);

  /** \fn DynamicalSystem* getDynamicalSystemOnNumber(int)
   *  \brief allows to get one specific DynamicalSystem, with the only identifier of this Dynamical System
   *  \param the identifier of the DynamicalSystem to get
   *  \return the DS having as identifier
   */
  DynamicalSystem* getDynamicalSystemOnNumber(int);

  /** \fn bool hasDynamicalSystemNumber(int)
   *  \brief determines if the DynamicalSystem which number is given in parameter exists
   *  \param the identifier of the DynamicalSystem to get
   *  \return bool
   */
  bool hasDynamicalSystemNumber(int);

  /** \fn inline int getInteractionVectorSize()
   *  \brief allows to get the size of the vector of Interactions
   *  \return int : the size of interactionVector
   */
  inline int getInteractionVectorSize() const
  {
    return this->interactionVector.size();
  }

  /** \fn vector<Interaction*> getInteractions(void)
   *  \brief allows to get all the Interaction of the NonSmoothDynamicalSystem problem
   *  \return the vector of Interaction
   */
  inline vector<Interaction*> getInteractions(void) const
  {
    return this->interactionVector;
  };

  /** \fn Interaction* getInteraction(int)
   *  \brief allows to get one specific Interaction, with the only identifier of this Interaction
   *  \param the identifier of the Interaction to get
   *  \return the Interaction having as identifier, the string given in parameter
   */
  Interaction* getInteraction(int);

  /** \fn Interaction* getInteractionOnNumber(int)
   *  \brief allows to get one specific Interaction, with his place in the interactionVector
   *  \param the identifier of the Interaction to get
   *  \return the Interaction at the right place in the interactionVector
   */
  Interaction* getInteractionOnNumber(int);

  /** \fn bool hasInteractionNumber(int)
   *  \brief determines if the Interaction which number is given in parameter exists
   *  \param the identifier of the Interaction to get
   *  \return bool
   */
  bool hasInteractionNumber(int);


  /** \fn void setBVP(bool)
   *  \brief allows to set the NonSmoothDynamicalSystem to BVP, else it is IVP
   *  \param bool : the value to set to BVP
   */
  inline void setBVP(const bool bvp)
  {
    this->BVP = bvp;
  };

  /** \fn void setBVP(bool)
   *  \brief allows to set the NonSmoothDynamicalSystem to IVP, else it is BVP
   *  \param bool : the value to set to IVP
   */
  inline void setIVP(const bool ivp)
  {
    this->BVP = !ivp;
  };

  /** \fn void setDynamicalSystems(vector<DynamicalSystem*>)
   *  \brief allows to set all the DynamicalSystems of the NonSmoothDynamicalSystem
   *  \param vector<DynamicalSystem> : the vector to set
   */
  inline void setDynamicalSystems(const vector<DynamicalSystem*> dsVect)
  {
    this->DSVector = dsVect;
  } ;

  /** \fn void setInteractions(vector<Interaction*>)
   *  \brief allows to set all the Interactions of the NonSmoothDynamicalSystem
   *  \param vector<Interaction> : the vector to set
   */
  inline void setInteractions(const vector<Interaction*> interVect)
  {
    this->interactionVector = interVect;
  };

  /** \fn inline NSDSXML* getNSDSXML()
   *  \brief allows to get the NSDSXML* of the NonSmoothDynamicalSystem problem
   *  \return the pointer on NSDSXML of the NonSmoothDynamicalSystem
   */
  inline NSDSXML* getNSDSXML() const
  {
    return this->nsdsxml;
  }

  /** \fn inline void setNSDSXML( NSDSXML *nsdsxml )
   *  \brief allows to set the NSDSXML* of the NonSmoothDynamicalSystem
   *  \param NSDSXML* : the pointer on NSDSXML* to set
   */
  inline void setNSDSXML(NSDSXML *nsdsxml)
  {
    this->nsdsxml = nsdsxml;
  }


  /** \fn vector<EqualityConstraint*> getEqualityConstraints(void)
   *  \brief allows to get all the EqualityConstraint of the NonSmoothDynamicalSystem problem
   *  \return the vector of EqualityConstraint
   */
  inline vector<EqualityConstraint*> getEqualityConstraints(void) const
  {
    return this->ecVector;
  };

  /** \fn EqualityConstraint* getEqualityConstraint(int)
   *  \brief allows to get one specific EqualityConstraint, with its place in the vector of EqualityConstraint
   *  \param int : the place of the EqualityConstraint in the vector of EqualityConstraint of the NonSmoothDynamicalSystem
   *  \return EqualityConstraint* : ecVector[ i ] EqualityConstraint
   */
  EqualityConstraint* getEqualityConstraint(int);

  /** \fn void setEqualityConstraints(vector<EqualityConstraint*>)
   *  \brief allows to set all the EqualityConstraints of the NonSmoothDynamicalSystem
   *  \param vector<EqualityConstraint*> : the vector to set
   */
  inline void setEqualityConstraints(const vector<EqualityConstraint*> ecVect)
  {
    this->ecVector = ecVect;
  };

  /** \fn void addEqualityConstraint(EqualityConstraint*)
   *  \brief allows to add the EqualityConstraint of the NonSmoothDynamicalSystem
   *  \param EqualityConstraint* : the EqualityConstraint to add
   */
  void addEqualityConstraint(EqualityConstraint*);

  /////////////////////////

  /** \fn void saveNSDSToXML()
   *  \brief copy the data of the NonSmoothDynamicalSystem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSDSToXML();

  /** \fn void display()
   *  \brief display the data of the Non Smooth Dynamical System
   */
  void display() const;

  /** \fn void addDynamicalSystem(DynamicalSystem*)
   *  \brief allows to add the DynamicalSystem to the NonSmoothDynamicalSystem
   *  \param DynamicalSystem* : the DynamicalSystem to add
   */
  void addDynamicalSystem(DynamicalSystem*);//, BoundaryCondition*);

  /** \fn void addInteraction(Interaction*)
   *  \brief allows to add the Interaction of the NonSmoothDynamicalSystem
   *  \param Interaction : the Interaction to add
   */
  void addInteraction(Interaction*);

  /** \fn Interaction* addInteraction(int number, int nInter, vector<int>* status, vector<DynamicalSystem*>*)
   *  \brief allow to add an Interaction to the NonSmoothDynamicalSystem
   *  \param int : the number of the Interaction
   *  \param int : the size of the y vector of the Interaction
   *  \param vector<int> : the status of the interaction
   */
  Interaction* addInteraction(int number = -1, int nInter = -1, vector<int>* status = NULL,
                              vector<DynamicalSystem*>* dsConcerned = NULL);

  /** \fn double nsdsConvergenceIndicator() const
   *  \brief calculate an indicator that gives convergence information for the DSs
   *  \return a double
   */
  double nsdsConvergenceIndicator() ;

protected:
  /** \fn void linkNSDSXML()
   *  \brief makes the links between the DSXMLs, InteractionXMLs of the NSDSXML of the NonSmoothDynamicalSystem and the DynamicalSystems, Interactions
   */
  void linkNSDSXML();

  /** \fn void fillNSDSWithNSDSXML()
   *  \brief uses the NSDSXML of the NonSmoothDynamicalSystem to fill the fields of this NonSmoothDynamicalSystem
   *  \exception RuntimeException
   */
  void fillNSDSWithNSDSXML();


private:
  /** TRUE if the NonSmoothDynamicalSystem is a boundary value problem*/
  bool BVP;

  /** contains Dynamic systems of the simulation */
  vector<DynamicalSystem*> DSVector;

  /** contains the Interactions */
  vector<Interaction*> interactionVector;

  /** contains the EqualityConstraints */
  vector<EqualityConstraint*> ecVector;

  /** the XML object linked to the NonSmoothDynamicalSystem to read XML data */
  NSDSXML *nsdsxml;
};

#endif
