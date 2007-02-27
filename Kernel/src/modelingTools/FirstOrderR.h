/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

/*! \file FirstOrderR.h
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.h"

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) Apr 27, 2004
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f[
 * y = h(X,t,\lambda,Z)\\
 * R = g(\lambda,t,Z)
 * \f]
 *  X, Z, R corresponds to DynamicalSystem variables.
 *  If DS1 and DS2 are involved in the linked Interaction, then X =[x1 x2], Z=[z1 z2] ...
 *  \f$ y \ and \ \lambda \f$ are specific variables of the interaction (see this class for more details).
 *  h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 *
 * h <=> output
 *
 * g <=> input
 *
 */
class FirstOrderR : public Relation
{

protected:

  /** type of the FirstOrderR */
  std::string  firstOrderType;

  /** FirstOrderR plug-in to compute y(x,t) - id="output".
   *  @param the size of the vector x.
   *  @param x : the pointer to the first element of the vector x.
   *  @param time : current time.
   *  @param the size of the vectors y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param[in,out]  y : the pointer to the first element of the vector y.
   *  @param the size of the vectors z.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  void (*computeOutputPtr)(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*);

  /** FirstOrderR plug-in to compute r(lambda,t) - id="input".
   *  @param sizeY : the size of the vector y and lambda.
   *  @param lambda : the pointer to the first element of the vector lambda.
   *  @param time : current time.
   *  @param[in,out] r : the pointer to the first element of the vector r.
   *  @param the size of the vectors z and lambda.
   *  @param[in,out] z : a vector of user-defined parameters.
   */
  void (*computeInputPtr)(unsigned int, const double*, double, double*, unsigned int, double*);

  /** default constructor
   *  \param a string that gives the type of the relation (optional)
   */
  FirstOrderR(const std::string& = "R");

public:

  /** xml constructor
   *  \param FirstOrderRXML* : the XML object corresponding
   *  \param a string that gives the type of the relation (optional)
   */
  FirstOrderR(RelationXML*, const std::string& = "R");

  /** data constructor
   *  \param a string with computeOutput function name.
   *  \param a string with computeInput function name.
   */
  FirstOrderR(const std::string&, const std::string&);

  /** destructor
   */
  virtual ~FirstOrderR();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  virtual void initialize();

  /** allows to get the type of the FirstOrderR
   *  \return string : the type of the FirstOrderR
   */
  inline const std::string  getFirstOrderRelationType() const
  {
    return firstOrderType;
  }

  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0);

  /** default function to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(double, unsigned int = 0);

  /** default function to compute r
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double, unsigned int);

  /** allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeOutputFunction(const std::string&, const std::string&);

  /** allow to set a specified function to compute output
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeInputFunction(const std::string&, const std::string&);

  /** main relation members display
   */
  virtual void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderR* convert(Relation *r);
};

#endif
