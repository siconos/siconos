/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/*! \file Relation.h
*/

#ifndef RELATION_H
#define RELATION_H

#include "SiconosConst.h"
#include "Interaction.h"
#include "RelationXML.h"
#include "DSInputOutput.h"
#include "check.h"

const std::string RELATION = "Relation";
const std::string LINEARTIRELATION = "LinearTIR";
const std::string LAGRANGIANRELATION = "LagrangianR";
const std::string LAGRANGIANLINEARRELATION = "LagrangianLinearR";

class Interaction;
class RelationXML;
class DSInputOutput;


//! General Non Linear Relation (Base class for Relations).
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
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

  /** Flags to know if pointers have been allocated inside constructors or not */
  std::map<std::string, bool> isAllocatedIn;

  /** Parameters list, last argument of plug-in functions. What are those parameters depends on user´s choice.
   *  This a list of pointer to SimpleVector. Each one is identified thanks to a key which is the plug-in name.
   * A flag is also added in the isAllocatedIn map to check inside-class memory allocation for this object.*/
  std::map<std::string, SimpleVector*> parametersList;

  /** computes y
   *  \param unsigned int sizeX : size of vector x
   *  \param double* x : the pointer to the first element of the vector x
   *  \param double* time : the current time
   *  \param unsigned int sizeY : size of vector y and lambda.
   *  \param double* lambda : the pointer to the first element of the vector lambda
   *  \param unsigned int sizeU : size of vector u
   *  \param double* u : the pointer to the first element of the vector u
   *  \param double* y : the pointer to the first element of the vector y (in-out parameter)
   *  \param double* param   : a vector of user-defined parameters
   */
  void (*computeOutputPtr)(const unsigned int, const double*, const double*, const unsigned int, const double*, const unsigned int, const double*, double*, double*);

  /** computes r
  *  \param unsigned int sizeY : size of vector y and lambda.
  *  \param double* lambda : the pointer to the first element of the vector lambda
  *  \param double* time : the current time
  *  \param double* r : the pointer to the first element of the vector r (in-out parameter)
  *  \param double* param   : a vector of user-defined parameters
  */
  void (*computeInputPtr)(const unsigned int, const double*, const double*, double*, double*);

  /** init parameter vector corresponding to id to a SimpleVector* of size 1
  *  \param a string, id of the plug-in
  */
  void initParameter(const std::string);

public:

  /** default constructor
  *  \param a string that gives the type of the relation (optional)
  */
  Relation(const std::string = "Relation");

  /** xml constructor
  *  \param RelationXML* : the XML object corresponding
  *  \param a string that gives the type of the relation (optional)
  */
  Relation(RelationXML*, const std::string = "Relation");

  /** copy constructor
  *  \param a relation to copy
  *  warning: the interaction link is not copied, set a new one!
  */
  Relation(const Relation&);

  /** destructor
  */
  virtual ~Relation();

  /** initialize the relation (check sizes, memory allocation ...)
  */
  virtual void initialize();

  /** allows to get the RelationXML* of the Relation
  *  \return a pointer on the RelationXML of the Relation
  */
  inline RelationXML* getRelationXML()
  {
    return relationxml;
  }

  /** allows to set the RelationXML* of the Relation
  *  \param RelationXML* : the pointer to set
  */
  inline void setRelationXML(RelationXML *rxml)
  {
    relationxml = rxml;
  }

  /** allows to get the Interaction which contains this Relation
  *  \return a pointer on an Interaction
  */
  inline Interaction* getInteractionPtr() const
  {
    return interaction;
  }

  /** set the Interaction which contains this Relation
  */
  inline void setInteractionPtr(Interaction* i)
  {
    interaction = i;
  }

  /** allows to get the type of the Relation
  *  \return string : the type of the Relation
  */
  inline const std::string  getType() const
  {
    return relationType;
  }

  /** get the name of computeInput function
  *  \return a string
  */
  inline const std::string getComputeInputName() const
  {
    return computeInputName;
  }

  /** get the name of computeOutput function
  *  \return a string
  */
  inline const std::string getComputeOutputName() const
  {
    return computeOutputName;
  }

  /** allows to get all the DSInputOutput of the Relation
  *  \return the vector of DSInputOutput
  */
  std::vector<DSInputOutput*> getDSInputOutputs(void);

  /** allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
  *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the Relation
  *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
  */
  DSInputOutput* getDSInputOutput(const unsigned int);

  /** allows to set all the DSInputOutputs of the Relation
  *  \param vector<DSInputOutput*> : the vector to set
  */
  void setDSInputOutputs(std::vector<DSInputOutput*>);

  /** allows to add the DSInputOutput to the Relation
  *  \param DSInputOutput* : the DSInputOutput to add
  */
  void addDSInputOutput(DSInputOutput*);

  // -- parametersList --

  /** get the full map of parameters
  *  \return a map<string,SimpleVector*>
  */
  inline std::map<std::string, SimpleVector*> getParameters() const
  {
    return parametersList;
  };

  /** get the vector of parameters corresponding to plug-in function named id
  *  \return a SimpleVector
  */
  inline const SimpleVector getParameter(const std::string id)
  {
    return *(parametersList[id]);
  };

  /** get the pointer to the vector of parameters corresponding to plug-in function named id
  *  \return a pointer on a SimpleVector
  */
  inline SimpleVector* getParameterPtr(const std::string id)
  {
    return parametersList[id];
  };

  /** set the map for parameters
  *  \param a map<string, SimpleVector*>
  */
  void setParameters(const std::map<std::string, SimpleVector*>&);

  /** set vector corresponding to plug-in function named id to newValue
  *  \param a SimpleVector
  *  \param a string
  */
  void setParameter(const SimpleVector&, const std::string);

  /** set vector corresponding to plug-in function named id to newPtr (!! pointer link !!)
  *  \param a pointer to SimpleVector
  *  \param a string
  */
  void setParameterPtr(SimpleVector *, const std::string);

  /** default function to compute y
  *  \param double : current time
  *  \param unsigned int: number of the derivative to compute, optional, default = 0.
  */
  virtual void computeOutput(const double, const unsigned int = 0);

  /** default function to compute y for the free state
  *  \param double : current time
  *  \param unsigned int: number of the derivative to compute, optional, default = 0.
  *  \exception RuntimeException
  */
  virtual void computeFreeOutput(const double, const unsigned int = 0);

  /** default function to compute r
  *  \param double : current time
  *  \param unsigned int: "derivative" order of lambda used to compute input
  */
  virtual void computeInput(const double, const unsigned int);

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeOutputFunction(const std::string, const std::string);

  /** allow to set a specified function to compute output
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeInputFunction(const std::string, const std::string);

  /** main relation members display
  */
  virtual void display() const;

  /** copy the data of the Relation to the XML tree
  */
  virtual void saveRelationToXML() const;

};

#endif // RELATION_H
