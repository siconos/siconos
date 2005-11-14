/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef RELATION_H
#define RELATION_H

#include "SiconosConst.h"
#include "Interaction.h"
#include "RelationXML.h"
#include "DSInputOutput.h"
#include "check.h"

class Interaction;
class RelationXML;
class DSInputOutput;

/** \class Relation
 *  \brief general non linear relation class.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 27, 2004
 *
 *    This class provides tools to define and describe relations of the type:
 * \f[
 * y = h(x,t,\lambda,u,...)
 * R = g(\lambda,t)
 * \f]
 *  x, u, R are DynamicalSystem variables.
 *  x being the dof vector for all the DS involved in the interaction that owns the current relation.
 *  u is a control term (see DynamicalSystem class), R the input due to non-smooth behavior.
 *   \f[ y and \lambda \f] are specific variables of the interaction (see this class for more details).
 * h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 * h <=> output
 * g <=> input
 *
 * They MUST be plugged -> to default plug-in functions if nothing else specified.
 */

class Relation
{

private:

  /** contains a link to the DSInputOutput of the DynamicalSystems */
  std::vector<DSInputOutput*> dsioVector;

protected:

  /** type of the Relation */
  std::string  relationType;

  /** the Interaction which contains this Relation */
  Interaction *interaction;

  /** the object linked this Relation to read XML data */
  RelationXML *relationxml;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* contains the name of the plugin used to compute g function */
  std::string  computeInputName;
  /* contains the name of the plugin used to compute h function */
  std::string  computeOutputName;

  /* Boolean variables to check if h and g are plugged to external functions or not
   *  - always true for general Relations
   *  - false by default for LinearTIR, but may be set to true by user, using setComputeOutput/Input functions
   *  - always false for Lagrangian ( "overloaded" with h(q,...) and G(q,...) )
   *  Note that these variables are only useful for derived classes.
   */
  bool isOutputPlugged;
  bool isInputPlugged;

  /** Parameters list, argument of plug-in. What are those parameters depends on user´s choice.
   *  The order corresponds to the one of the plug-in list below :
   *  computeOutput, computeInput
   */
  std::vector<SimpleVector*> parametersList0; // -> Size = 2
  std::deque<bool> isParametersList0AllocatedIn;

  /** \fn void (*computeOutputPtr)(const unsigned int* sizeX, const double* x, const double* time,
                                   const unsigned int* sizeY, const double* lambda,
           const unsigned int* sizeU, const double* u, double* y, double* param);
   *  \brief computes y
   *  \param unsigned int* sizeX : size of vector x
   *  \param double* x : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param unsigned int* sizeY : size of vector y and lambda.
   *  \param double* lambda : the pointer to the first element of the vector lambda
   *  \param unsigned int* sizeU : size of vector u
   *  \param double* u : the pointer to the first element of the vector u
   *  \param double* y : the pointer to the first element of the vector y (in-out parameter)
   *  \param double* param   : a vector of user-defined parameters
   */
  void (*computeOutputPtr)(const unsigned int*, const double*, const double*, const unsigned int*, const double*, const unsigned int*, const double*, double*, double*);

  /** \fn void (*computeInputPtr)(const unsigned int* sizeY, const double* lambda, const double* time, double* r, double* param);
   *  \brief computes r
   *  \param unsigned int* sizeY : size of vector y and lambda.
   *  \param double* lambda : the pointer to the first element of the vector lambda
   *  \param double* time : the current time
   *  \param double* r : the pointer to the first element of the vector r (in-out parameter)
   *  \param double* param   : a vector of user-defined parameters
   */
  void (*computeInputPtr)(const unsigned int*, const double*, const double*, double*, double*);

public:

  /** \fn Relation(Interaction* =NULL)
   *  \brief default constructor
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  Relation(Interaction* = NULL);

  /** \fn Relation(RelationXML*)
   *  \brief xml constructor
   *  \param RelationXML* : the XML object corresponding
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  Relation(RelationXML*, Interaction* = NULL);

  /** \fn Relation(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  Relation(const Relation&, Interaction* = NULL);

  /** \fn ~Relation()
   *  \brief destructor
   */
  virtual ~Relation();

  /** \fn inline RelationXML* getRelationXML()
   *  \brief allows to get the RelationXML* of the Relation
   *  \return a pointer on the RelationXML of the Relation
   */
  inline RelationXML* getRelationXML()
  {
    return relationxml;
  }

  /** \fn inline void setRelationXML(RelationXML *rxml)
   *  \brief allows to set the RelationXML* of the Relation
   *  \param RelationXML* : the pointer to set
   */
  inline void setRelationXML(RelationXML *rxml)
  {
    relationxml = rxml;
  }

  /** \fn Interaction* getInteractionPtr()
   *  \brief allows to get the Interaction which contains this Relation
   *  \return a pointer on an Interaction
   */
  inline Interaction* getInteractionPtr() const
  {
    return interaction;
  }

  /** \fn void setInteractionPtr(Interaction* i)
   *  \brief set the Interaction which contains this Relation
   */
  virtual inline void setInteractionPtr(Interaction* i)
  {
    interaction = i;
  }

  /** \fn inline string getType()
   *  \brief allows to get the type of the Relation
   *  \return string : the type of the Relation
   */
  inline const std::string  getType() const
  {
    return relationType;
  }

  /** \fn inline string getComputeInputName()
   *  \brief get the name of computeInput function
   *  \return a string
   */
  inline const std::string  getComputeInputName() const
  {
    return computeInputName;
  }

  /** \fn inline string getComputeOutputName()
   *  \brief get the name of computeOutput function
   *  \return a string
   */
  inline const std::string  getComputeOutputName() const
  {
    return computeOutputName;
  }

  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the Relation
   *  \return the vector of DSInputOutput
   */
  std::vector<DSInputOutput*> getDSInputOutputs(void);

  /** \fn DSInputOutput* getDSInputOutput(const int&)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the Relation
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the Relation
   *  \param vector<DSInputOutput*> : the vector to set
   */
  void setDSInputOutputs(std::vector<DSInputOutput*>);

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the Relation
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput*);

  // -- parametersList --

  /** \fn vector<SimpleVector*> getParametersListVector(unsigned int & index) const
   *  \brief get the parameter list at position index
   *  \return SimpleVector
   */
  virtual inline std::vector<SimpleVector*> getParametersListVector() const
  {
    return parametersList0;
  }

  /** \fn  const SimpleVector getParametersList(const unsigned int & index) const
   *  \brief get the parameter list at position index
   *  \return SimpleVector
   */
  virtual inline const SimpleVector getParametersList(const unsigned int & index) const
  {
    return *(parametersList0[index]);
  }

  /** \fn SimpleVector* getParametersListPtr(const unsigned int & index) const
   *  \brief get the pointer on the parameter list at position index
   *  \return pointer on a SimpleVector
   */
  virtual inline SimpleVector* getParametersListPtr(const unsigned int & index) const
  {
    return parametersList0[index];
  }

  /** \fn void setParametersListVector(const std::vector<SimpleVector*>& newVector)
   *  \brief set vector parameterList0 to newVector
   *  \param vector<SimpleVector*>
   */
  virtual void setParametersListVector(const std::vector<SimpleVector*>&);

  /** \fn void setParametersList (const SimpleVector& newValue, const unsigned int & index)
   *  \brief set the value of parameterList0[index] to newValue
   *  \param SimpleVector newValue
   */
  virtual void setParametersList(const SimpleVector&, const unsigned int &);

  /** \fn void setParametersListPtr(SimpleVector* newPtr, const unsigned int & index)
   *  \brief set parametersList0[index] to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  virtual  void setParametersListPtr(SimpleVector *newPtr, const unsigned int & index);

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeOutput(const double&);

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeFreeOutput(const double&);

  /** \fn void computeInput(double time);
   *  \brief default function to compute r
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeInput(const double&);

  /** \fn void setComputeOutputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeOutputFunction(const std::string&, const std::string&);

  /** \fn void setComputeInputFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeInputFunction(const std::string&, const std::string&);

  /** \fn  void display() const
   * \brief main relation members display
   */
  virtual void display() const;

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  virtual void saveRelationToXML() const;


};

#endif // RELATION_H
