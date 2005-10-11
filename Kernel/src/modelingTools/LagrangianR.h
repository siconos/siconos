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
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.h"
#include "LagrangianRXML.h"
#include "LagrangianDS.h"

/** \class LagrangianR
 *  \brief Lagrangian (Non Linear) Relation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date Apr 27, 2004
 *  This class provides tools to describe non linear relation of the type:
 *
 * \todo : write complete equations for NL Lagrangian relations
 *
 *  using plug-in mechanism for function H definition.
 *
 * Recall:
 *   - output is vector Y, ie y and all its derivatives from 1 to r-1, r being the relative degree (see Interaction)
 *   - input is vector Lambda (see Interaction)
 *
 * Some rules:
 *   - h must be a plug-in
 *   - jacobianQH can be either a fixed matrix or a plug-in, to fill this matrix in.
 *
 * Remark: at the time, we consider that h is a function of q (and not velocity).
 * That means jacobianVelocity is not yet implemented.
 *
 */
class LagrangianRXML;

class LagrangianR : public Relation
{

protected:

  // --- h ---

  /* bool to check whether h is a plug-in or not */
  /* Warning: for LagrangianR, h is always a plug-in but not necessary for LagrangianLinearR derived class */
  /* Recall: h(q) = Hq + b in derived class */
  bool isHPlugin;
  /* Name of the plugin used to compute h */
  std::string  hFunctionName;

  /** \fn void (*computeHPtr)(const unsigned int* sizeOfq, const double* time,
   *                          const unsigned int* sizeOfy, const double* q, const double* y);
   * \brief computes y = h(q,t)
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: time
   * \param unsigned int: size of output vector y
   * \param double*: vector q (composite of all DS q-vectors)
   * \param double*: output vector y (in-out parameter)
   */
  void (*computeHPtr)(const unsigned int*, const double*, const unsigned int*, const double*, double*);

  // --- jacobianH ---

  /* Jacobian Matrix of function h and a boolean to know whether it has been allocated or not */
  /* at the time, we on;y use jacobian compare to q, index 0 of stl vectors */
  std::vector<SiconosMatrix*> jacobianH;
  std::vector<bool> isJacobianHAllocatedIn;
  // Remark: stl vector usefull for future implementation with jacobian compare to velocity and so on

  /* bool to check whether jacobian is a plug-in or not. See remark above concerning the use of stl vector */
  std::vector<bool> isJacobianHPlugin;

  /* Name of the plugin used to compute jacobian of H compare to q */
  /* See remark above concerning the use of stl vector */
  std::vector<std::string>  jacobianHFunctionName;

  /** \fn void (*computeJacobianQHPtr)(const unsigned int* sizeOfq, const double* time,
   *                                  const unsigned int* sizeOfy, const double* q, const double* jacob);
   * \brief computes jacobian compare to q of function h
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: time
   * \param unsigned int: size of output vector y
   * \param double*: vector q (composite of all DS q-vectors)
   * \param double*: output vector jacobianH (in-out parameter)
   */
  void (*computeJacobianQHPtr)(const unsigned int*, const double*, const unsigned int*, const double*, double*);

  /** class for plug-in management (open, close librairy...) */
  SiconosSharedLibrary cShared;

public:

  /** \fn LagrangianR(Interaction* =NULL);
   * \brief default constructor
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianR(Interaction* = NULL);

  /** \fn void LagrangianR(RelationXML*, Interaction* =NULL)
   *  \brief constructor from xml file
   *  \param relationXML
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
   */
  LagrangianR(RelationXML*, Interaction* = NULL);

  /** \fn void LagrangianR(const string& computeInput,const string& computeOutput, Interaction* =NULL)
   *  \brief constructor from a set of data
   *  \param string : the name of the plugin for computeInput
   *  \param string : the name of the plugin for computeOutput
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
   */
  LagrangianR(const std::string&, const std::string&, Interaction* = NULL);

  /** \fn LagrangianR(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  LagrangianR(const Relation &, Interaction* = NULL);

  ~LagrangianR();

  // -- JacobianQH --

  /** \fn  const SiconosMatrix getJacobianQH() const
   *  \brief get the value of JacobianQH
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianQH() const
  {
    return *(jacobianH[0]);
  }

  /** \fn SiconosMatrix* getJacobianQHPtr() const
   *  \brief get JacobianQH
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQHPtr() const
  {
    return jacobianH[0];
  }

  /** \fn void setJacobianQH (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianQH to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianQH(const SiconosMatrix&);

  /** \fn void setJacobianQHPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianQH to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianQHPtr(SiconosMatrix *newPtr);


  /** \fn void setComputeHFunction(const string pluginPath, const string functionName&)
   *  \brief to set a specified function to compute function h(q)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeHFunction(const std::string &, const std::string &);

  /** \fn void setComputeJacobianHFunction(const string pluginPath, const string functionName&)
   *  \brief to set a specified function to compute jacobian of h compare to q
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQHFunction(const std::string &, const std::string &);


  /** \fn  std::string getHFunctionName() const
   *  \brief get name of function that computes h
   *  \return a string
   */
  inline const std::string getHFunctionName() const
  {
    return hFunctionName;
  }

  /** \fn  std::vector<string> getJacobianHFunctionName() const
   *  \brief get names of functions that compute jacobian
   *  \return a vector of strings
   */
  inline const std::vector<std::string> getJacobianHFunctionName() const
  {
    return jacobianHFunctionName;
  }

  /** \fn const bool isHPlugged() const
   *  \brief bool that shows if h is a plug-in or not
   *  \return a bool
   */
  inline const bool isHPlugged() const
  {
    return isHPlugin;
  };

  /** \fn const vector<bool> isJacobianHPlugged() const
  *  \brief list of bool that shows if jacobianH is a plug-in or not
  *  \return a vector of bool
  */
  inline const std::vector<bool> isJacobianHPlugged() const
  {
    return isJacobianHPlugin;
  };

  /** \fn void computeH(const double & time);
   * \brief to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(const double &);

  /** \fn void computeJacobianQH(const double & time);
   * \brief to compute Jacobian of h(q,v,t) compare to q using plug-in mechanism
   * \param: double, current time
   */
  void computeJacobianQH(const double &);

  /** \fn void computeOutput(double time);
   *  \brief to compute output
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeOutput(const double&);

  /** \fn void computeFreeOutput(double time);
   *  \brief to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeFreeOutput(const double&);

  /** \fn void computeInput(double time);
   *  \brief to compute p
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeInput(const double&);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  void saveRelationToXML();

  /** \fn LagrangianR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianR* convert(Relation *r);

  /** \fn  void display() const
   * \brief main relation members display
   */
  virtual void display() const;

};

#endif // LAGRANGIANRELATION_H
