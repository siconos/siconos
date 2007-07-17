/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
/*! \file LagrangianR.h

*/
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.h"

class Relation;
class DynamicalSystem;
class LagrangianRXML;

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr2)(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr3)(unsigned int, const double*, unsigned int, double*, unsigned int, double*);

/**Pointer to function - Plug-in utilities*/
typedef void (*FPtr4)(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*);

/** Lagrangian (Non Linear) Relation (generic interface)
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date Apr 27, 2004
 *
 * Relations for Lagrangian Dynamical Systems. This class is only an interface for specific (Linear, Scleronomous ...)
 * Lagrangian Relations (see derived classes).
 *
 *  Class name = type+subType.
 *
 * If y = h(...), all the gradients of are handled by G object.
 * For example, G[0] = \f$ \nabla_q h(q,...) \f$.
 *
 * In corresponding derived classes, h and Gi are connected to plug-in functions (user-defined).
 *
 */
class LagrangianR : public Relation
{

protected:

  /** To define the type of constraints (scleronomic ...), ie the variables on which depend h and G*/
  std::string LagrangianRelationType;

  /** G matrices of gradients of h. */
  VectorOfMatrices G;

  /** default constructor
   */
  LagrangianR(const std::string& = "Lagrangian");

  /** constructor from xml file
   *  \param relationXML
   *  \param string: relation subType
   */
  LagrangianR(RelationXML*, const std::string&);

public:

  /** destructor
   */
  virtual ~LagrangianR();

  /** xml utility: used to read G[i] in xml file
      \param LagrangianRXML, the xml pointer
      \param i, index of required G
  */
  void readGInXML(LagrangianRXML *, unsigned int);

  /** initialize the relation (check sizes, memory allocation ...)
   */
  virtual void initialize();

  /** get the type of constraints of the relation (scleronomic ...)
   *  \return a string
   */
  inline const std::string getLagrangianRelationType() const
  {
    return LagrangianRelationType;
  }

  /** initialize G matrices or components specific to derived classes.
   */
  virtual void initComponents();

  // -- G --

  /** get the vector of matrices G
   *  \return vector<SiconosMatrix*>
   */
  inline VectorOfMatrices getGVector() const
  {
    return G;
  }

  /** set the values of G[i]: no copy but links with pointers!!
   *  \param a VectorOfMatrices.
   */
  void setGVector(const VectorOfMatrices&);

  /** get matrix G[index]
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getG(unsigned int  index = 0) const
  {
    return *(G[index]);
  }

  /** get a pointer on matrix G[index]
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getGPtr(unsigned int index = 0) const
  {
    return G[index];
  }

  /** set the value of G[index] to newValue (copy)
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in G vector
   */
  void setG(const SiconosMatrix&, unsigned int = 0);

  /** set G[index] to pointer newPtr (pointer link)
   *  \param SiconosMatrix * newPtr
   *  \param unsigned int: index position in G vector
   */
  void setGPtr(SiconosMatrix *newPtr, unsigned int = 0);

  /** to set a specified function to compute function h(q,...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  virtual void setComputeHFunction(const std::string& , const std::string&);

  /** to set a specified function to compute G(q, ...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \param unsigned int: the index of G that must be computed (see introduction of this class for details on indexes)
   */
  virtual void setComputeGFunction(const std::string& , const std::string&  , unsigned int  = 0);

  /** to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  virtual void computeH(double);

  /** to compute G using plug-in mechanism. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  virtual void computeG(double, unsigned int = 0);

  /** to compute output
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0) = 0;

  /** to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(double, unsigned int = 0) = 0;

  /** to compute p
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double, unsigned int) = 0;

  /** copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** main relation members display
   */
  virtual void display() const;

};

#endif // LAGRANGIANRELATION_H
