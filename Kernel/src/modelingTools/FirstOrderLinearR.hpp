/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
/*! \file FirstOrderLinearR.hpp

 */
#ifndef FirstOrderLinearR_H
#define FirstOrderLinearR_H

#include "FirstOrderR.hpp"
/** Pointer to function used for plug-in for matrix-type operators (C,F etc) */
typedef void (*FOMatPtr1)(double, unsigned int, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for square matrix-type operators (D) */
typedef void (*FOMatPtr2)(double, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for vector-type operators (e) */
typedef void (*FOVecPtr)(double, unsigned int, double*, unsigned int, double*);

/** First Order Linear Relation

\author SICONOS Development Team - copyright INRIA
\version 3.0.0.
\date Apr 15, 2007

Linear Relation for First Order Dynamical Systems:

\f{eqnarray}
y &=& C(t,z)x(t) + F(t,z)z + D(t,z)\lambda + e(t,z) \\

R &=& B(t,z) \lambda
\f}

All coefficients can be plugged or not. Use isPlugged[name] to check if name ( = "C", "F" etc) is plugged.

Note: the connections (pointers equalities) between C, D, B and jacobianH and jacobianG of FirstOrderR class are done during initialize.

 */
class FirstOrderLinearR : public FirstOrderR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearR);

  SP::SiconosMatrix _F;
  SP::SiconosVector _e;

public:

  /** default constructor, protected
  */
  FirstOrderLinearR();

  /** Constructor with C and B plugin names
  \param Cname the plugin name for computing the C matrix
  \param Bname the plugin name for computing the B matrix
  **/
  FirstOrderLinearR(const std::string& Cname, const std::string& Bname);

  /** Constructor with all plugin names
  \param Cname the plugin name for computing the C matrix
  \param Dname the plugin name for computing the D matrix
  \param Fname the plugin name for computing the F matrix
  \param Ename the plugin name for computing the e vector
  \param Bname the plugin name for computing the B matrix
  **/
  FirstOrderLinearR(const std::string& Cname, const std::string& Dname, const std::string& Fname, const std::string& Ename, const std::string& Bname);

  /** create the Relation from a set of data
  *  \param C the C matrix
  *  \param B the B matrix
  */
  FirstOrderLinearR(SP::SiconosMatrix C, SP::SiconosMatrix B);

  /** create the Relation from a set of data
  *  \param C the C matrix
  *  \param D the D matrix
  *  \param F the F matrix
  *  \param e the e matrix
  *  \param B the B matrix
  */
  FirstOrderLinearR(SP::SiconosMatrix C, SP::SiconosMatrix D, SP::SiconosMatrix F, SP::SiconosVector e, SP::SiconosMatrix B);

  /** destructor
  */
  ~FirstOrderLinearR() {};

  // GETTERS/SETTERS

  /** set a specified function to compute the matrix C
  *  \param pluginPath the complete path to the plugin
  *  \param functionName the function name to use in this plugin
  */
  void setComputeCFunction(const std::string& pluginPath, const std::string& functionName)
  {
    setComputeJachxFunction(pluginPath,  functionName);
  }

  /** set a specified function to compute the matrix D
  *  \param pluginPath the complete path to the plugin
  *  \param functionName the function name to use in this plugin
  */
  void setComputeDFunction(const std::string& pluginPath, const std::string& functionName)
  {
    setComputeJachlambdaFunction(pluginPath,  functionName);
  }

  /** set a specified function to compute the matrix B
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   */
  void setComputebFunction(const std::string& pluginPath, const std::string& functionName)
  {
    setComputeJacglambdaFunction(pluginPath,  functionName);
  }

  /** set D to pointer newPtr
  *  \param newPtr the D matrix
  */
  inline void setDPtr(SP::SiconosMatrix newPtr)
  {
    _jachlambda = newPtr;
  }

  /** set F to pointer newPtr
  *  \param newPtr the F matrix
  */
  inline void setFPtr(SP::SiconosMatrix newPtr)
  {
    _F = newPtr;
  }

  /** get F
  *  \return F matrix
  */
  inline SP::SiconosMatrix F() const
  {
    return _F;
  }

  /** get e
  *  \return e matrix
  */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /** Function to compute matrix C
  */
  virtual void computeC(double time, Interaction& inter);

  /** Function to compute matrix F
  */
  virtual void computeF(double time, Interaction& inter);

  /** Function to compute matrix D
  */
  virtual void computeD(double time, Interaction& inter);

  /** Function to compute vector e
  */
  virtual void computeE(double time, Interaction& inter);

  /** Function to compute matrix B
  */
  virtual void computeb(double time, Interaction& inter);

  /** default function to compute h
  *  \param double : current time
  */
  virtual void computeh(double time, Interaction& inter);

  /** default function to compute g
  *  \param double : current time
  */
  virtual void computeg(double time, Interaction& inter);

  /** default function to compute y
  *  \param double: not used
  *  \param unsigned int: not used
  */
  virtual void computeOutput(double time, Interaction& inter, unsigned int = 0);

  /** default function to compute r
  *  \param double : not used
  *  \param unsigned int: not used
  */
  void computeInput(double time, Interaction& inter, unsigned int = 0);

  /** initialize the relation (check sizes, memory allocation ...)
  *  \param inter the interaction that owns this relation
  */
  virtual void initialize(Interaction & inter);

  /** print the data to the screen
  */
  void display() const;

  /** determine if the Relation is linear
  * \return true if the relation is linear.
  */
  virtual bool isLinear()
  {
    return true;
  }
  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearR)

#endif
