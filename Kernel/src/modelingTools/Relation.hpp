/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file Relation.hpp
  \brief General interface for relations.
*/

#ifndef RELATION_H
#define RELATION_H

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "RuntimeException.hpp"
#include "Tools.hpp"
#include "SiconosPointers.hpp"
#include "PluginTypes.hpp"
#include "RelationNamespace.hpp"
#include "PluggedObject.hpp"
#include "DynamicalSystemsSet.hpp"
//#include "Interaction.hpp"

class SiconosVector;
class SimpleMatrix;
class Interaction;

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
 * - computeInput(...) to compute non-smooth DS part (r or p) using
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
 *  - a pointer to an xml object
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
  /** Plug-in to compute \f$ \nabla_{\lambda} h(..)\f$
   */
  SP::PluggedObject _pluginJachlambda;

  /** Plug-in to compute g(...)
   */
  SP::PluggedObject _pluging;

  /** Plug-in to compute \f$ \nabla_\lambda g(lambda,t,z)\f$
   */
  SP::PluggedObject _pluginJacLg;
  /** Plug-in to compute f.$
   */
  SP::PluggedObject _pluginf;
  /** Plug-in to compute e.$
   */
  SP::PluggedObject _plugine;
  /** To initialize all the plugin functions with NULL.
   */
  virtual void zeroPlugin();

  SP::SiconosMatrix _jachlambda;

  /** type of the Relation: FirstOrder or Lagrangian */
  RELATION::TYPES _relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  RELATION::SUBTYPES _subType;

  /** the object linked this Relation to read XML data */
  SP::RelationXML _relationxml;

  /** basic constructor
   *  \param a string that gives the type of the relation
   *  \param a string that gives the subtype of the relation
   */
  Relation(RELATION::TYPES, RELATION::SUBTYPES);

  /** xml constructor
   *  \param RelationXML* : the XML object corresponding
   *  \param a string that gives the type of the relation
   *  \param a string that gives the subtype of the relation
   */
  Relation(SP::RelationXML, RELATION::TYPES, RELATION::SUBTYPES);


private:

  /** default constructor => private, no copy nor pass-by-value
   */
  Relation();

  /** copy constructor => private, no copy nor pass-by-value.
   */
  Relation(const Relation&);

  /** Assignment  => private, forbidden
   */
  Relation& operator=(const Relation&);

public:


  /** destructor
   */
  virtual ~Relation();

  /** To get the RelationXML* of the Relation
   *  \return a pointer on the RelationXML of the Relation
   */
  inline SP::RelationXML getRelationXML()
  {
    return _relationxml;
  }

  /** To set the RelationXML* of the Relation
   *  \param RelationXML* : the pointer to set
   */
  inline void setRelationXML(SP::RelationXML rxml)
  {
    _relationxml = rxml;
  }

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

  /** To get the name of h plugin
   *  \return a string
   */
  const std::string gethName() const ;

  /** To get the name of g plugin
   *  \return a string
   */
  const std::string getgName() const;

  /** To get the name of Jach[i] plugin
   *  \return a string
   */
  virtual const std::string getJachxName() const ;

  /** To get the name of Jacg[i] plugin
   *  \return a string
   */
  virtual const std::string getJacgName(unsigned int) const;

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
  virtual void setComputeEFunction(const std::string& pluginPath, const std::string& functionName);

  /** To set a plug-in function to compute the jacobian according to x of the input
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  virtual void setComputeJacglambdaFunction(const std::string& pluginPath, const std::string& functionName);

  /** initialize the relation (check sizes, memory allocation ...)
      \param inter the interaction using this relation
  */
  virtual void initialize(Interaction& inter) = 0;

  /** default function to compute h
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeh(const double time, Interaction& inter) = 0;

  /** default function to compute g
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeg(const double time, Interaction& inter);

  /** default function to compute jacobianG according to lambda
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacglambda(const double time, Interaction& inter) = 0;

  /** compute all the H Jacobian
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJach(const double time, Interaction& inter) = 0;

  /* compute all the G Jacobian
   *  \param time the current time
   *  \param inter the interaction using this relation
   */
  virtual void computeJacg(const double time, Interaction& inter) = 0;


  /** default function to compute y
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param derivativeNumber number of the derivative to compute (optional, default = 0)
   */
  virtual void computeOutput(const double time, Interaction& inter, unsigned int derivativeNumber = 0) = 0;

  /** default function to compute r
   *  \param time the current time
   *  \param inter the interaction using this relation
   *  \param level the input "derivative" order of lambda used to compute input
   */
  virtual void computeInput(const double time, Interaction& inter, unsigned int level = 0) = 0;

  virtual inline SP::SiconosMatrix jachlambda() const
  {
    return _jachlambda;
  }


  virtual SP::SiconosMatrix C() const = 0;


  /**
   * return true if the relation is linear.
   */

  virtual bool isLinear()
  {
    return false;
  }

  /** main relation members display
   */
  virtual void display() const;

  /** copy the data of the Relation to the XML tree
   */
  virtual void saveRelationToXML() const;

  /** Check if _pluginh is correctly set */
  inline bool ishPlugged() const
  {
    return _pluginh->isPlugged();
  };

  /** Check if _pluginJachx is correctly set */
  inline bool isJachxPlugged() const
  {
    return _pluginJachx->isPlugged();
  };

  /** Check if _pluginJachlambda is correctly set */
  inline bool isJachlambdaPlugged() const
  {
    return _pluginJachlambda->isPlugged();
  };

  /** Check if _pluging is correctly set */
  inline bool isgPlugged() const
  {
    return _pluging->isPlugged();
  };

  /** Check if _pluginJacLg is correctly set */
  inline bool isJacLgPlugged() const
  {
    return _pluginJacLg->isPlugged();
  };

  /** Check if _pluginf is correctly set */
  inline bool isfPlugged() const
  {
    return _pluginf->isPlugged();
  };

  /** Check if _plugine is correctly set */
  inline bool isePlugged() const
  {
    return _plugine->isPlugged();
  };

  /** Get _pluginh */
  inline SP::PluggedObject getPluginh() const
  {
    return _pluginh;
  };

  /** Get _pluginJachx */
  inline SP::PluggedObject getPluginJachx() const
  {
    return _pluginJachx;
  };

  /** Get _pluginJachlambda */
  inline SP::PluggedObject getPluginJachlambda() const
  {
    return _pluginJachlambda;
  };

  /** Get _pluging */
  inline SP::PluggedObject getPluging() const
  {
    return _pluging;
  };

  /** Get _pluginJacLg */
  inline SP::PluggedObject getPluginJacLg() const
  {
    return _pluginJacLg;
  };

  /** Get _pluginf */
  inline SP::PluggedObject getPluginf() const
  {
    return _pluginf;
  };

  /** Get _plugine */
  inline SP::PluggedObject getPlugine() const
  {
    return _plugine;
  };
  /** visitors hook
   */
  virtual void preparNewtonIteration(Interaction& inter)
  {
    ;
  };
  VIRTUAL_ACCEPT_VISITORS(Relation);

};

#endif // RELATION_H
