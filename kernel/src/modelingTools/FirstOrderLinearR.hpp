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

The following operators can be plugged: \f$ B(t,z), C(t,z), D(t,z), e(t,z), F(t,z)\f$

 */
class FirstOrderLinearR : public FirstOrderR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearR);

  /** initialize the relation (check sizes, memory allocation in workV and workM ...)
  *  \param inter Interaction using this Relation
  *  \param DSlink
  *  \param workV
  *  \param workM
  */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

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
  FirstOrderLinearR(SP::SimpleMatrix C, SP::SimpleMatrix B);

  /** create the Relation from a set of data
  *  \param C the C matrix
  *  \param D the D matrix
  *  \param F the F matrix
  *  \param e the e matrix
  *  \param B the B matrix
  */
  FirstOrderLinearR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e, SP::SimpleMatrix B);

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
  void setComputeBFunction(const std::string& pluginPath, const std::string& functionName)
  {
    setComputeJacglambdaFunction(pluginPath,  functionName);
  }

  /** Function to compute the matrix C
   * \param time the current time
   * \param z the auxiliary input vector
   * \param C the C matrix
  */
  void computeC(double time, SiconosVector& z, SimpleMatrix& C);

  /** Function to compute the matrix D
   * \param time the current time
   * \param z the auxiliary input vector
   * \param D the D matrix
  */
  void computeD(double time, SiconosVector& z, SimpleMatrix& D);

  /** Function to compute the matrix F
   * \param time the current time
   * \param z the auxiliary input vector
   * \param F the F matrix
  */
  void computeF(double time, SiconosVector& z, SimpleMatrix& F);

  /** Function to compute the vector e
   * \param time the current time
   * \param z the auxiliary input vector
   * \param e the e vector
  */
  void computee(double time, SiconosVector& z, SiconosVector& e);

  /** Function to compute the matrix B
   * \param time the current time
   * \param z the auxiliary input vector
   * \param B the B matrix
  */
  void computeB(double time, SiconosVector& z, SimpleMatrix& B);

  /** default function to compute h
  *  \param time current time
  *  \param workV
  *  \param workM
  *  \param x XXX
  *  \param z XXX
  *  \param lambda
  *  \param y value of h
  */
  void computeh(double time, VectorOfVectors& workV, VectorOfSMatrices& workM,
                BlockVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y);

  /** default function to compute g
  *  \param time current time
  *  \param workM
  *  \param lambda XXX
  *  \param z XXX
  *  \param r non-smooth input
  */
  void computeg(double time, VectorOfSMatrices& workM, SiconosVector& lambda, SiconosVector& z, BlockVector& r);

  /** default function to compute y
  *  \param time current time
  *  \param inter Interaction using this Relation
  *  \param interProp
  *  \param level not used
  */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** default function to compute r
  *  \param time current time
  *  \param inter Interaction using this Relation
  *  \param interProp
  *  \param level not used
  */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** print the data to the screen
  */
  void display() const;

  /** set e
  *  \param  newe the new value of e
  */
  inline void setePtr(SP::SiconosVector newe)
  {
    _e = newe;
  }

  /** get e
  *  \return e matrix
  */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /** determine if the Relation is linear
  * \return true if the relation is linear.
  */
  virtual bool isLinear()
  {
    return true;
  }

  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp) {};
  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp) {};


  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearR)

#endif
