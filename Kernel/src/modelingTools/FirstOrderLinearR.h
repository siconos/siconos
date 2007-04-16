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
/*! \file FirstOrderLinearR.h

*/
#ifndef FirstOrderLinearR_H
#define FirstOrderLinearR_H

#include "FirstOrderR.h"

/** Pointer to function used for plug-in for matrix-type operators (C,F etc) */
typedef void (*FOMatPtr1)(double, unsigned int, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for square matrix-type operators (D) */
typedef void (*FOMatPtr2)(double, unsigned int, double*, unsigned int, double*);

/** Pointer to function used for plug-in for vector-type operators (e) */
typedef void (*FOVecPtr)(double, unsigned int, double*, unsigned int, double*);


/** Linear Time Invariant Relation, derived from class FirstOrderR
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date Apr 27, 2004
 *
 *  Linear Relation for First Order Dynamical Systems:
 *
 * \f[
 * y = C(t,z)x(t) + F(t,z)z + D(t,z)\lambda + e(t,z) \\
 *
 * R = B(t,z) \lambda
 * \f]
 *
 * All coefficients can be plugged or not. Use isPlugged[name] to check if name ( = "C", "F" etc) is plugged.
 *
 * Note: the connections (pointers equalities) between C, D, B and jacobianH and jacobianG of FirstOrderR class are done during initialize.
 *
 */
class FirstOrderLinearR : public FirstOrderR
{

protected:

  /** C */
  SiconosMatrix* C;
  /** D*/
  SiconosMatrix* D;
  /** F*/
  SiconosMatrix* F;
  /** e*/
  SiconosVector* e;
  /**  B*/
  SiconosMatrix* B;

  /** Plug-in list to compute C(t,z)
   * @param time: current time
   * @param rowOfC, number of rows of in-out matrix C
   * @param colOfC, number of columns of in-out matrix C
   * @param[in,out] C : pointer to the first element of C
   * @param[in] sizeOfZ: size of vector z
   * @param[in,out] z: a vector of user-defined parameters
   */
  FOMatPtr1 CPtr;

  /** Plug-in to compute D(t,z)
   * @param time: current time
   * @param rowOfD, number of rows of in-out square matrix D
   * @param[in,out] D : pointer to the first element of D
   * @param[in] sizeOfZ: size of vector z
   * @param[in,out] z: a vector of user-defined parameters
   */
  FOMatPtr2 DPtr;

  /** Plug-in to compute F(t,z)
   * @param time: current time
   * @param rowOfF, number of rows of in-out matrix F
   * @param[in,out] F : pointer to the first element of F
   * @param[in] sizeOfZ: size of vector z
   * @param[in,out] z: a vector of user-defined parameters
   */
  FOMatPtr2 FPtr;

  /** Plug-in to compute e(t,z)
   * @param time: current time
   * @param sizeOfE, size of in-out vector e
   * @param[in,out] e : pointer to the first element of e
   * @param[in] sizeOfZ: size of vector z
   * @param[in,out] z: a vector of user-defined parameters
   */
  FOVecPtr ePtr;

  /** Plug-in to compute B(t,z)
   * @param time: current time
   * @param rowOfB, number of rows of in-out matrix B
   * @param colOfB, number of columns of in-out matrix B
   * @param[in,out] B : pointer to the first element of B
   * @param[in] sizeOfZ: size of vector z
   * @param[in,out] z: a vector of user-defined parameters
   */
  FOMatPtr1 BPtr;

  /** protected function used to initilized isPlugged map
      \param : a bool, value for all flags.
  */
  void initPluginFlags(bool);

  /** protected function used to initilized isAllocated map
      \param : a bool, value for all flags.
  */
  void initAllocationFlags(bool);

  /** Default (private) constructor
   */
  FirstOrderLinearR();

public:

  /** xml constructor
   *  \param FirstOrderLinearRXML* : the XML object corresponding
   */
  FirstOrderLinearR(RelationXML*);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   *  \param SiconosMatrix : the matrix B
   *  \exception RuntimeException
   */
  FirstOrderLinearR(const SiconosMatrix& , const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SimpleVector  : e
   *  \param SiconosMatrix : B
   *  \exception RuntimeException
   */
  FirstOrderLinearR(const SiconosMatrix& , const SiconosMatrix& ,
                    const SiconosMatrix& , const SimpleVector& ,
                    const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param pointer to SiconosMatrix : the matrix C
   *  \param pointer to SiconosMatrix : the matrix B
   *  \exception RuntimeException
   */
  FirstOrderLinearR(SiconosMatrix* , SiconosMatrix*);

  /** create the Relation from a set of data
   *  \param pointer to SiconosMatrix : C
   *  \param pointer to SiconosMatrix : D
   *  \param pointer to SiconosMatrix : F
   *  \param pointer to SimpleVector  : e
   *  \param pointer to SiconosMatrix : B
   *  \exception RuntimeException
   */
  FirstOrderLinearR(SiconosMatrix* , SiconosMatrix* ,
                    SiconosMatrix* , SimpleVector* ,
                    SiconosMatrix*);

  /** destructor
   */
  virtual ~FirstOrderLinearR();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  void initialize();

  // GETTERS/SETTERS

  // -- C --

  /** get the value of C
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getC() const
  {
    return *C;
  }

  /** get C
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getCPtr() const
  {
    return C;
  }

  /** set the value of C to newValue
   *  \param SiconosMatrix newValue
   */
  void setC(const SiconosMatrix&);

  /** set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setCPtr(SiconosMatrix *);

  /** set a specified function to compute the matrix C
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeCFunction(const std::string& , const std::string&);

  // -- D --

  /** get the value of D
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getD() const
  {
    return *D;
  }

  /** get D
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getDPtr() const
  {
    return D;
  }

  /** set the value of D to newValue
   *  \param SiconosMatrix newValue
   */
  void setD(const SiconosMatrix&);

  /** set D to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setDPtr(SiconosMatrix *);

  /** set a specified function to compute the matrix D
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeDFunction(const std::string& , const std::string&);

  // -- F --

  /** get the value of F
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getF() const
  {
    return *F;
  }

  /** get F
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getFPtr() const
  {
    return F;
  }

  /** set the value of F to newValue
   *  \param SiconosMatrix newValue
   */
  void setF(const SiconosMatrix&);

  /** set F to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setFPtr(SiconosMatrix *) ;

  /** set a specified function to compute the matrix F
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeFFunction(const std::string& , const std::string&);

  // -- e --

  /** get the value of e
   *  \return SimpleVector
   */
  inline const SimpleVector getE() const
  {
    return *e;
  }

  /** get e
   *  \return pointer on a SimpleVector
   */
  inline SiconosVector* getEPtr() const
  {
    return e;
  }

  /** set the value of e to newValue
   *  \param SiconosVector newValue
   */
  void setE(const SiconosVector&);

  /** set E to pointer newPtr
   *  \param  SiconosVector* newPtr
   */
  void setEPtr(SiconosVector*);

  /** set a specified function to compute the vector e
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeEFunction(const std::string& , const std::string&);

  // -- B --

  /** get the value of B
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getB() const
  {
    return *B;
  }

  /** get B
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getBPtr() const
  {
    return B;
  }

  /** set the value of B to newValue
   *  \param SiconosMatrix newValue
   */
  void setB(const SiconosMatrix&);

  /** set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setBPtr(SiconosMatrix *) ;

  /** set a specified function to compute the matrix B
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   */
  void setComputeBFunction(const std::string& , const std::string&);

  /** Function to compute matrix C
   */
  void computeC(double);

  /** Function to compute matrix F
   */
  void computeF(double);

  /** Function to compute matrix D
   */
  void computeD(double);

  /** Function to compute vector e
   */
  void computeE(double);

  /** Function to compute matrix B
   */
  void computeB(double);

  /** Computes y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(double, unsigned int = 0);

  /** Computes yFree AND save it into y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(double = 0, unsigned int = 0);

  /** Computes lambda
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(double, unsigned int);

  /** default function to compute jacobianH (from FirstOrderR)
   *  \param double : current time
   *  \param index for jacobian (0: jacobian according to x, 1 according to lambda)
   */
  void computeJacobianH(double, unsigned int);

  /** default function to compute jacobianG according to lambda (from FirstOrderR)
   *  \param double : current time
   *  \param index for jacobian: at the time only one possible jacobian => i = 0 is the default value .
   */
  void computeJacobianG(double, unsigned int = 0);

  /** copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** print the data to the screen
   */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderLinearR* convert(Relation *r);
};

#endif
