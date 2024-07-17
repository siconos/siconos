/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <memory>
#include <string>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosSerialization.hpp"  // For ACCEPT_SERIALIZATION
#include "SiconosVector.hpp"

namespace siconos::algebra {
// class SiconosMatrix;
// class SiconosVector;
// class BlockVector;
}  // namespace siconos::algebra

namespace siconos::plugins {
class PluggedObject;
}

namespace siconos::internal {

struct SiconosVisitor;
}

namespace siconos::modeling {

class Interaction;

/** List of possible Relations types*/
enum class RelationType {
  /** First Order */
  FirstOrder,
  /** Lagrangian */
  Lagrangian,
  /** Lagrangian */
  NewtonEuler
};

/** List of possible Relations subtypes*/
enum class RelationSubType {
  /** non linear */
  NonLinearR,
  /** linear */
  LinearR,
  /** Linear and time invariant */
  LinearTIR,
  /** Scleronomous (lagrangian only) */
  ScleronomousR,
  /** Rheonomous (lagrangian only) */
  RheonomousR,
  /** Compliant (lagrangian only) */
  CompliantR,
  /** CompliantLinearTIR (lagrangian only) */
  CompliantLinearTIR,
  /** */
  Type1R,
  /** */
  Type2R
};

/**
   General Non Linear Relation (Abstract Base class for Relations).

   The present class is an interface to all relations and provides
   tools to define and describe them.

   A relation is a link between global variables of the Dynamical
   Systems and some local ones, named y and lambda; belonging to one
   and only one Interaction.

   All relations are specified by their type (First order or Lagrangian)
   accessed by getType() and their sub-type (linear, scleronomous ...),
   returned by getSubType().

   A relation provides functions to compute:

   - a function computeOutput() that updates y using dynamical systems global
   variables,
   - a function computeInput() that updates non-smooth dynamical systems parts
   (e.g. r or p) using  \f$ \lambda \f$ .

*/
class Relation {
 protected:
  ACCEPT_SERIALIZATION(Relation);

  /** Plug-in to compute h(...)
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginh{nullptr};

  /** Plug-in to compute \f$ \nabla_x h(..) \f$
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJachx{nullptr};
  /** Plug-in to compute \f$ \nabla_z h(..) \f$
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJachz{nullptr};
  /** Plug-in to compute \f$ \nabla_{\lambda} h(..) \f$
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJachlambda{nullptr};

  /** Plug-in to compute g(...)
   */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluging{nullptr};

  /** Plug-in to compute \f$ \nabla_\lambda g \f$  */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJacglambda{nullptr};

  /** Plug-in to compute \f$ \nabla_x g \f$  */
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginJacgx{nullptr};

  /** Plug-in to compute f*/
  std::shared_ptr<siconos::plugins::PluggedObject> _pluginf{nullptr};

  /** Plug-in to compute e*/
  std::shared_ptr<siconos::plugins::PluggedObject> _plugine{nullptr};

  /** To initialize all the plugin functions with nullptr.
   */
  virtual void _zeroPlugin();

  /** type of the Relation: FirstOrder or Lagrangian */
  RelationType _relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  RelationSubType _subType;

  /** base and only constructor
   *
   *  \param type type of the relation
   *  \param subtype subtype of the relation
   */
  Relation(RelationType type, RelationSubType subtype)
      : _relationType(type), _subType(subtype) {
    _zeroPlugin();
  };

 private:
  // Rule of five ...
  Relation() = delete;
  Relation(const Relation &) = delete;
  Relation(Relation &&) = delete;
  Relation &operator=(const Relation &) = delete;
  Relation &operator=(Relation &&) = delete;

 public:
  /** destructor */
  virtual ~Relation() noexcept = default;

  /** \return the type of the Relation (FirstOrder or Lagrangian) */
  inline RelationType getType() const { return _relationType; }

  /** \return the subType of the Relation */
  inline RelationSubType getSubType() const { return _subType; }

  /** To set a plug-in function to compute output function h
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputehFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute   \f$  \nabla_x h(..) \f$
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeJachxFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute  \f$ \nabla_z h(..) \f$
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeJachzFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute  \f$ \nabla_{\lambda} h(..) \f$
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeJachlambdaFunction(const std::string &pluginPath,
                                    const std::string &functionName);

  /** To set a plug-in function to compute input function g
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputegFunction(const std::string &pluginPath, const std::string &functionName);
  /** To set a plug-in function to compute input function F
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeFFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute input function E
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeEFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute the jacobian of  \f$ g \f$  w.r.t. x
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeJacgxFunction(const std::string &pluginPath, const std::string &functionName);

  /** To set a plug-in function to compute the jacobian of  \f$ g \f$  w.r.t. \f$ \lambda \f$
   *
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputeJacglambdaFunction(const std::string &pluginPath,
                                    const std::string &functionName);

  /** initialize the relation (check sizes, memory allocation ...)
   *
   *  \param inter the interaction using this relation
   */
  virtual void initialize(Interaction &inter) = 0;

  /** check sizes of the relation specific operators.
   *
   *  \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction &inter) = 0;

  /** compute all the H Jacobian
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(double time, Interaction &inter) = 0;

  /** compute all the G Jacobian
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param interProp
   */
  virtual void computeJacg(double time, Interaction &inter) = 0;

  /** default function to compute y
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber number of the derivative to compute (optional,
   *  default = 0)
   */
  virtual void computeOutput(double time, Interaction &inter,
                             unsigned int derivativeNumber = 0) = 0;
  /** default function to compute r
   *
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param level the input "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double time, Interaction &inter, unsigned int level = 0) = 0;

  virtual std::shared_ptr<siconos::algebra::SiconosMatrix> C() const = 0;

  virtual std::shared_ptr<siconos::algebra::SiconosMatrix> H() const = 0;

  /**
     return true if the relation is linear.

     \return bool
   */
  virtual bool isLinear() { return false; }

  /**
      return true if the relation requires the computation of residu

      \return true if residu are required, false otherwise
   */
  virtual bool requireResidu() { return false; }

  /** main relation members display
   */
  virtual void display() const;

  /**
      Get _pluginh

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginh() const {
    return _pluginh;
  };

  /**
      Get _pluginJachx

      \return a shared pointer to the plugin
   */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginJachx() const {
    return _pluginJachx;
  };

  /**
      Get _pluginJachlambda

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginJachlambda() const {
    return _pluginJachlambda;
  };

  /**
      Get _pluging

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluging() const {
    return _pluging;
  };

  /**
      Get _pluginJacglambda

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginJacLg() const {
    return _pluginJacglambda;
  };

  /**
      Get _pluginf

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPluginf() const {
    return _pluginf;
  };

  /**
      Get _plugine

      \return a shared pointer to the plugin
  */
  inline std::shared_ptr<siconos::plugins::PluggedObject> getPlugine() const {
    return _plugine;
  };

  virtual void accept(std::shared_ptr<siconos::internal::SiconosVisitor>) const = 0;

  // VIRTUAL_ACCEPT_VISITORS(Relation);
};
}  // namespace siconos::modeling

#endif  // RELATION_H
