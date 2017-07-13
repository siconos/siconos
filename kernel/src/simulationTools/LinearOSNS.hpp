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
/*! \file LinearOSNS.hpp
  \brief Linear Complementarity Problem formulation and solving
*/

#ifndef LinearOSNS_H
#define LinearOSNS_H

#include "OneStepNSProblem.hpp"

/** stl vector of double */
typedef std::vector<double> MuStorage;
TYPEDEF_SPTR(MuStorage)

/** Base (abstract) class for linear non-smooth problems

    \author SICONOS Development Team - copyright INRIA
    \date (Creation) November 13, 2010

    Base class for linear non-smooth problems, usually in the form:

    \f[
    w =  q + M z
    \f]
    where
    - \f$ w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknowns,
    - \f$ M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$

    examples: LCP, FrictionContact ...

*/
class LinearOSNS : public OneStepNSProblem
{
protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LinearOSNS);

  /** vector w of a LinearOSNS system */
  SP::SiconosVector _w;

  /** vector z of a LinearOSNS system */
  SP::SiconosVector _z;

  /** matrix M of a LinearOSNS system */
  SP::OSNSMatrix _M;

  /** vector q of a LinearOSNS system */
  SP::SiconosVector _q;

  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage
      (embedded into OSNSMatrix) */
  int _numericsMatrixStorageType;

  /** a boolean to decide if _w and _z vectors are initialized with
      previous values of Y and Lambda when a change occurs in problem
      size */
  bool _keepLambdaAndYState;

  /** nslaw effects : visitors experimentation
   */
  struct _TimeSteppingNSLEffect;
  struct _EventDrivenNSLEffect;
  struct _NSLEffectOnSim;
  friend struct _TimeSteppingNSLEffect;
  friend struct _EventDrivenNSLEffect;
  friend struct _NSLEffectOnSim;

  /** default constructor (private)
   */
  LinearOSNS() ;

public:

  /** constructor from data
      \param numericsSolverId the numerics_solver identifier
  */
  LinearOSNS(const int numericsSolverId);

  /** destructor
   */
  virtual ~LinearOSNS() {};

  // --- W ---
  /** copy of the current value of vector w
      \return a SiconosVector
   */
  inline const SiconosVector getW() const
  {
    return *_w;
  }

  /** current w vector (pointer link)
      \return pointer on a SiconosVector
  */
  inline SP::SiconosVector w() const
  {
    return _w;
  }

  /** set w vector (pointer link)
      \param newPtr the new SP::SiconosVector
  */
  inline void setWPtr(SP::SiconosVector newPtr)
  {
    _w = newPtr;
  }

  // --- Z ---
  /** copy of the current value of vector z
      \return a SiconosVector
   */
  inline const SiconosVector getz() const
  {
    return *_z;
  }

  /** current z vector (pointer link)
      \return pointer on a SiconosVector
  */
  inline SP::SiconosVector z() const
  {
    return _z;
  }

  /** set z vector (pointer link)
      \param newPtr the new SP::SiconosVector
  */
  inline void setzPtr(SP::SiconosVector newPtr)
  {
    _z = newPtr;
  }

  /** M matrix (pointer link)
      \return pointer on a OSNSMatrix
   */
  inline SP::OSNSMatrix M() const
  {
    return _M;
  }

  /** set M to pointer newPtr
      \param newM the new M matrix
   */
  inline void setMPtr(SP::OSNSMatrix newM)
  {
    _M = newM;
  }

  /** get the value of q, the constant vector in the LinearOSNS
      \return SiconosVector
   */
  inline const SiconosVector getQ() const
  {
    return *_q;
  }

  /** get q, the the constant vector in the LinearOSNS
      \return pointer on a SiconosVector
  */
  inline SP::SiconosVector q() const
  {
    return _q;
  }

  /** set q to pointer newPtr
      \param newQ the new q vector
   */
  inline void setQPtr(SP::SiconosVector newQ)
  {
    _q = newQ;
  }

  /** get the type of storage used for M
      \return int (0: dense, 1:sparse)
   */
  inline int getMStorageType() const
  {
    return _numericsMatrixStorageType;
  };

  /** set which type of storage will be used for M
   * \warning this function does not allocate any memory for M,
   * it just sets an indicator for future use
   * \param i (0:dense, 1:sparse)
   */
  inline void setMStorageType(int i)
  {
    _numericsMatrixStorageType = i;
  };

  /** Memory allocation or resizing for z,w,q */
  void initVectorsMemory();

  /** initialize the _M matrix */
  virtual void initOSNSMatrix();
 
  /** To initialize the LinearOSNS problem(computes topology ...)
      \param sim the simulation owning this OSNSPB
  */
  virtual void initialize(SP::Simulation sim);

  /** compute extra-diagonal interactionBlock-matrix
   *  \param ed an edge descriptor
   */
  virtual void computeInteractionBlock(const InteractionsGraph::EDescriptor& ed);

  /** compute diagonal Interaction block
   * \param vd a vertex descriptor
   */
  virtual void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd);

  /** To compute a part of the "q" vector of the OSNS
      \param vertex, vertex (interaction) which corresponds to the considered block
      \param pos the position of the first element of yOut to be set
  */
  virtual void computeqBlock(InteractionsGraph::VDescriptor& vertex, unsigned int pos);

  /** compute vector q
   *  \param time the current time
   */
  virtual void computeq(double time);

  /** build problem coefficients (if required)
      \param time the current time
      \return true if succeeded
   */
  virtual bool preCompute(double time);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *  \param time the current time
   *  \return information about the solver convergence.
   */
  virtual int compute(double time) = 0;

  /** update interactions variables (y and lambda) according to current problem found solutions.
   */
  virtual void postCompute();

  /** print the data to the screen
   */
  virtual void display() const;

  /** choose initialisation behavior for w and z.
      \param val true: init w and z with previous values
      of y and lambda saved in interactions, false: init to 0.
  */
  void setKeepLambdaAndYState(bool val)
  {
    _keepLambdaAndYState = val ;
  }

  /* visitors hook */
  ACCEPT_STD_VISITORS();


};

#endif // LinearOSNS_H
