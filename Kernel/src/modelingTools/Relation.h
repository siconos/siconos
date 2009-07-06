/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
  \brief General interface for relations.
*/

#ifndef RELATION_H
#define RELATION_H

#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "Interaction.h"
#include "RuntimeException.h"
#include "Tools.hpp"
#include "Plugin.hpp"
#include "SiconosPointers.hpp"
#include "PluginTypes.hpp"
#include "RelationNamespace.hpp"
#include <boost/weak_ptr.hpp>

class SimpleVector;
class SimpleMatrix;

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

  /** type of the Relation: FirstOrder or Lagrangian */
  RELATION::TYPES relationType;

  /** sub-type of the Relation (exple: LinearTIR or ScleronomousR ...) */
  RELATION::SUBTYPES subType;

  /** The Interaction linked to this Relation */
  boost::weak_ptr<Interaction> interaction;

  /** True if h is plugged to a user-defined function */
  bool hPlugged;

  /** True if g is plugged to a user-defined function */
  bool gPlugged;

  /** Name of the plugin function used to compute h*/
  std::string hName;

  /** Name of the plugin function used to compute g*/
  std::string gName;

  /** A map of vectors, used to save links (pointers) to DS objects of
      the interaction */
  std::vector<SP::SiconosVector> data;

  /** the object linked this Relation to read XML data */
  SP::RelationXML relationxml;

  /** work vector for x */
  SP::SimpleVector workR;

  /** work vector for x */
  SP::SimpleVector workX;

  /** work vector for z */
  SP::SimpleVector workZ;

  /** work vector for y */
  SP::SimpleVector workY;

  /** work vector for lambda */
  SP::SimpleVector workL;

  /** The residu y of the newton iterations*/
  SP::SiconosVector mResiduy;
  /*value of h at the current newton iteration*/
  SP::SiconosVector mH_alpha;

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

  /** To get the pointer to the Interaction linked to the present
   *  Relation
   *  \return a pointer to Interaction.
   */
  inline SP::Interaction getInteractionPtr()
  {
    assert(!interaction.expired());
    return interaction.lock();
  }

  /** To set the pointer to the Interaction linked to the present
   *  Relation
   *  \param a pointer to Interaction.
   */
  inline void setInteractionPtr(SP::Interaction newInter)
  {
    interaction = boost::weak_ptr<Interaction>(newInter);
  }

  /** To get the RelationXML* of the Relation
   *  \return a pointer on the RelationXML of the Relation
   */
  inline SP::RelationXML getRelationXML()
  {
    return relationxml;
  }

  /** To set the RelationXML* of the Relation
   *  \param RelationXML* : the pointer to set
   */
  inline void setRelationXML(SP::RelationXML rxml)
  {
    relationxml = rxml;
  }

  /** To get the type of the Relation (FirstOrder or Lagrangian)
   *  \return the type of the Relation
   */
  inline const RELATION::TYPES  getType() const
  {
    return relationType;
  }

  /** To get the subType of the Relation
   *  \return the sub-type of the Relation
   */
  inline const RELATION::SUBTYPES  getSubType() const
  {
    return subType;
  }

  /** To get the name of h plugin
   *  \return a string
   */
  inline const std::string getHName() const
  {
    return hName;
  }

  /** To get the name of g plugin
   *  \return a string
   */
  inline const std::string getGName() const
  {
    return gName;
  }

  /** To get the name of JacH[i] plugin
   *  \return a string
   */
  virtual const std::string getJacHName(unsigned int) const
  {
    return "unamed";
  }

  /** To get the name of JacG[i] plugin
   *  \return a string
   */
  virtual const std::string getJacGName(unsigned int) const
  {
    return "unamed";
  }

  /** true if h is plugged
   *  \return a bool
   */
  inline const bool isHPlugged() const
  {
    return hPlugged;
  }

  /** true if g is plugged
   *  \return a bool
   */
  inline const bool isGPlugged() const
  {
    return gPlugged;
  }

  /** true if JacH[i] is plugged
   *  \return a bool
   */
  virtual const bool isJacHPlugged(unsigned int) const
  {
    return false;
  }

  /** true if JacG[i] is plugged
   *  \return a bool
   */
  virtual const bool isJacGPlugged(unsigned int) const
  {
    return false;
  }

  /** get matrix JacH[index]
   *  \return a SimpleMatrix
   */
  virtual const SimpleMatrix getJacH(unsigned int  index = 0) const = 0;

  /** get a pointer on matrix JacH[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual SP::SiconosMatrix getJacHPtr(unsigned int index = 0) const = 0;

  /** get matrix JacG[index]
   *  \return a SimpleMatrix
   */
  virtual const SimpleMatrix getJacG(unsigned int  index = 0) const
  {
    RuntimeException::selfThrow("Relation::getJacG() - not implemented for this type of relation (probably Lagrangian): " + getType());
    return SimpleMatrix(0, 0);
  }

  /** get a pointer on matrix JacG[index]
   *  \return a pointer on a SiconosMatrix
   */
  virtual SP::SiconosMatrix getJacGPtr(unsigned int index = 0) const
  {
    RuntimeException::selfThrow("Relation::getJacGPtr() - not implemented for this type of relation (probably Lagrangian): " + getType());
    return SP::SiconosMatrix();
  }

  /** Gets the number of computed jacobians for h
      \return an unsigned int.
  */
  virtual unsigned int getNumberOfJacobiansForH() const
  {
    return 0;
  }

  /** Gets the number of computed jacobian for g
      \return an unsigned int.
  */
  virtual unsigned int getNumberOfJacobiansForG() const
  {
    return 0;
  }

  /** To set a plug-in function to compute output function h
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  virtual void setComputeHFunction(const std::string&, const std::string&) = 0;

  /** To set a plug-in function to compute jacobianH
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index for jacobian (0: jacobian according to x, 1
   *  according to lambda)
   */
  virtual void setComputeJacobianHFunction(const std::string&, const std::string&, unsigned int = 0) = 0;

  /** To set a plug-in function to compute input function g
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  virtual void setComputeGFunction(const std::string&, const std::string&)
  {
    RuntimeException::selfThrow("Relation::setComputeGFunction() - not implemented for this type of relation (probably Lagrangian): " + getType());
  }

  /** To set a plug-in function to compute the jacobian according to x of the input
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void setComputeJacobianGFunction(const std::string&, const std::string&, unsigned int = 0)
  {
    RuntimeException::selfThrow("Relation::setComputeJacobianGFunction() - not implemented for this type of relation (probably Lagrangian): " + getType());
  }

  /** initialize the relation (check sizes, memory allocation ...)
      \param SP to Interaction: the interaction that owns this relation
  */
  virtual void initialize(SP::Interaction) = 0;
  virtual void initializeMemory() ;

  /** default function to compute h
   *  \param double : current time
   */
  virtual void computeH(double) = 0;

  /** default function to compute g
   *  \param double : current time
   */
  virtual void computeG(double t);

  /** default function to compute jacobianH
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  virtual void computeJacH(double, unsigned int) = 0;

  /** default function to compute jacobianG according to lambda
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible
   *  jacobian => i = 0 is the default value .
   */
  virtual void computeJacG(double, unsigned int)
  {
    ;//RuntimeException::selfThrow("Relation::computeJacG() - not implemented for this type of relation (probably Lagrangian): "+getType());
  }
  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute,
   *  optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0) = 0;

  /** default function to compute r
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to
   *  compute input
   */
  virtual void computeInput(double, unsigned int = 0) = 0;


  /**
   * return a SP on the C matrix.
   * The matrix C in the linear case, else it returns Jacobian of the output with respect to x.
   *
   */
  virtual SP::SiconosMatrix getCPtr();
  /**
   * return a SP on the D matrix.
   * The matrix D in the linear case, else it returns Jacobian of the output with respect to lambda.
   */
  virtual SP::SiconosMatrix getDPtr();
  /**
   * return a SP on the B matrix.
   * The matrix B in the linear case, else it returns Jacobian of the input with respect to lambda.
   */
  virtual SP::SiconosMatrix getBPtr();
  // --- Residu y functions

  inline const SP::SiconosVector getResiduY()
  {
    return mResiduy;
  }
  /*
   * Compute the residuY.
   *
   *
   */
  virtual void computeResiduY(double t);

  /*
   * Return H_alpha
   *
   */
  virtual const SP::SiconosVector getHalphaPtr()
  {
    return mH_alpha;
  };

  /*
   * Do the computation needed by for a Iteration.
   *
   */
  virtual void preparNewtonIteration()
  {
    ;
  };

  /**
   * return true if the relation is linear.
   */

  virtual bool isLinear()
  {
    return true;
  }

  /** main relation members display
   */
  virtual void display() const;

  /** copy the data of the Relation to the XML tree
   */
  virtual void saveRelationToXML() const;


  /** visitor hook for shared pointers
   */
  virtual void accept(SP::SiconosVisitor)
  {
    RuntimeException::selfThrow
    ("this class derived from Relation does not accept a visitor for shared pointers");
  };

  /** visitor hook for standard pointers
   */
  virtual void accept(SiconosVisitor&)
  {
    RuntimeException::selfThrow
    ("this class derived from Relation does not accept a visitor");
  }

};
#endif // RELATION_H
