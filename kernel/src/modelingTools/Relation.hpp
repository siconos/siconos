/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file Relation.hpp
  \brief General interface for relations.
*/

#ifndef RELATION_H
#define RELATION_H
#include "SiconosFwd.hpp"
#include "RelationNamespace.hpp"
#include "SiconosVisitor.hpp"
#include "SiconosAlgebraTypeDef.hpp"

/** General Non Linear Relation (Virtual Base class for Relations).
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *  A relation is a link between global variables of the Dynamical
 * Systems and some local ones, named y and lambda; belonging to one
 * and only one Interaction.
 *
 * The present class is an interface to all relations and provides
 * tools to define and describe them.
 *
 * Each relation must have the two following functions:
 *
 *  - computeOutput(...) to compute y using DS global variables.
 *  - computeInput(...) to compute non-smooth DS part (r or p) using
 *   lambda.
 *
 * Depending on the DS class and the link type, various relations (ie
 * derived classes) are available:
 *   - FirstOrder, for FirstOrderDS and derived classes.
 *   - Lagrangian, for LagrangianDS and derived classes.
 *
 *  The specific type (Linear, Scleronomous ...) is then given by the
 *  "subType". \n
 *
 * The relation holds also:
 *  - a VectorMap to handle links to DS variables (no copy!!!). Filled
 *    in during initialize.
 *
 */
class Relation
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Relation);

  /** Plug-in to compute h(...)
  */
  SP::PluggedObject _pluginh;

  /** Plug-in to compute \f$ \nabla_x h(..)\f$
   */
  SP::PluggedObject _pluginJachx;
  /** Plug-in to compute \f$ \nabla_z h(..)\f$
   */
  SP::PluggedObject _pluginJachz;
  /** Plug-in to compute \f$ \nabla_{\lambda} h(..)\f$
   */
  SP::PluggedObject _pluginJachlambda;

  /** Plug-in to compute g(...)
   */
  SP::PluggedObject _pluging;

  /** Plug-in to compute \f$ \nabla_\lambda g\f$ */
  SP::PluggedObject _pluginJacglambda;

  /** Plug-in to compute \f$ \nabla_x g\f$ */
  SP::PluggedObject _pluginJacgx;

  /** Plug-in to compute f*/
  SP::PluggedObject _pluginf;

  /** Plug-in to compute e*/
  SP::PluggedObject _plugine;
  /** To initialize all the plugin functions with NULL.
   */
  virtual void _zeroPlugin();

  /** type of the Relation: FirstOrder or Lagrangian */
  RELATION::TYPES _relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  RELATION::SUBTYPES _subType;

  /** basic constructor
   *  \param type type of the relation
   *  \param subtype subtype of the relation
   */
  Relation(RELATION::TYPES type, RELATION::SUBTYPES subtype);

private:

  /** default constructor => private, no copy nor pass-by-value
   */
  Relation();

  /** copy constructor => private, no copy nor pass-by-value.
   */
  Relation(const Relation&);

  /* Assignment  => private, forbidden
   */
  Relation& operator=(const Relation&);

public:


  /** destructor
   */
  virtual ~Relation(){};

  /** To get the type of the Relation (FirstOrder or Lagrangian)
   *  \return the type of the Relation
   */
  inline RELATION::TYPES  getType() const
  {
    return _relationType;
  }

  /** To get the subType of the Relation
   *  \return the sub-type of the Relation
   */
  inline RELATION::SUBTYPES  getSubType() const
  {
    return _subType;
  }

  /** To set a plug-in function to compute output function h
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputehFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute  \f$ \nabla_x h(..)\f$
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJachxFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute  \f$ \nabla_z h(..)\f$
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJachzFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute  \f$ \nabla_{\lambda} h(..)\f$
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJachlambdaFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute input function g
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputegFunction(const std::string& pluginPath, const std::string& functionName);
  /** To set a plug-in function to compute input function F
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeFFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute input function E
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeEFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute the jacobian of \f$g\f$ w.r.t. x
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJacgxFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute the jacobian of \f$g\f$ w.r.t. \f$\lambda\f$.
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJacglambdaFunction(const std::string& pluginPath, const std::string& functionName);

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   * \param DSlink the container of the link to DynamicalSystem attributes
   * \param workV work vectors
   * \param workM work matrices
   */
  virtual void initialize(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM) = 0;
  
  /** compute all the H Jacobian
   * \param time the current time
   * \param inter the interaction using this relation
   * \param interProp
   */
  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp) = 0;

  /* compute all the G Jacobian
   * \param time the current time
   * \param inter the interaction using this relation
   * \param interProp
   */
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp) = 0;


  /** default function to compute y
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param interProp
   *  \param derivativeNumber number of the derivative to compute (optional, default = 0)
   */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int derivativeNumber = 0) = 0;

  /** default function to compute r
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param interProp
   *  \param level the input "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0) = 0;

  virtual SP::SimpleMatrix C() const = 0;

  /** return true if the relation is linear.
      \return bool
   */
  virtual bool isLinear()
  {
    return false;
  }

  /** return true if the relation requires the computation of residu
      \return true if residu are required, false otherwise
   */
  virtual bool requireResidu()
  {
    return false;
  }

  /** main relation members display 
   */
  virtual void display() const;

  /** Get _pluginh 
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPluginh() const
  {
    return _pluginh;
  };

  /** Get _pluginJachx
      \return a shared pointer to the plugin
   */
  inline SP::PluggedObject getPluginJachx() const
  {
    return _pluginJachx;
  };

  /** Get _pluginJachlambda
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPluginJachlambda() const
  {
    return _pluginJachlambda;
  };

  /** Get _pluging
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPluging() const
  {
    return _pluging;
  };

  /** Get _pluginJacglambda
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPluginJacLg() const
  {
    return _pluginJacglambda;
  };

  /** Get _pluginf
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPluginf() const
  {
    return _pluginf;
  };

  /** Get _plugine
      \return a shared pointer to the plugin
  */
  inline SP::PluggedObject getPlugine() const
  {
    return _plugine;
  };
  /** visitors hook
   *  \param inter  interaction 
   *  \param interProp
   */
  virtual void prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp)
  {
    ;
  };
  VIRTUAL_ACCEPT_VISITORS(Relation);

};

#endif // RELATION_H
