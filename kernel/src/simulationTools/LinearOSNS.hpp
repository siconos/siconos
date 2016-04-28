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
    \version 3.0.0.
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
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LinearOSNS);


  /** contains the vector w of a LinearOSNS system */
  SP::SiconosVector _w;

  /** contains the vector z of a LinearOSNS system */
  SP::SiconosVector _z;

  /** contains the matrix M of a LinearOSNS system */
  SP::OSNSMatrix _M;

  /** contains the vector q of a LinearOSNS system */
  SP::SiconosVector _q;

  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage
      (embedded into OSNSMatrix) */
  int _MStorageType;

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
  /** get the value of w, the initial state of the DynamicalSystem
   *  \return a SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an
   *  lvalue => return SiconosVector
   */
  inline const SiconosVector getW() const
  {
    return *_w;
  }

  /** get w, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector w() const
  {
    return _w;
  }

  /** set the value of w to newValue
   *  \param newValue the new SiconosVector
   */
  void setW(const SiconosVector& newValue);

  /** set w to pointer newPtr
   *  \param newPtr the new SP::SiconosVector
   */
  inline void setWPtr(SP::SiconosVector newPtr)
  {
    _w = newPtr;
  }

  // --- Z ---
  /** get the value of z, the initial state of the DynamicalSystem
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an
   *  lvalue => return SiconosVector
   */
  inline const SiconosVector getz() const
  {
    return *_z;
  }

  /** get z, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector z() const
  {
    return _z;
  }

  /** set the value of z to newValue
   *  \param newValue a SiconosVector
   */
  void setz(const SiconosVector& newValue);

  /** set z to pointer newPtr
   *  \param newPtr a SP::SiconosVector
   */
  inline void setzPtr(SP::SiconosVector newPtr)
  {
    _z = newPtr;
  }

  // --- M ---

  /** get M
   *  \return pointer on a OSNSMatrix
   */
  inline SP::OSNSMatrix M() const
  {
    return _M;
  }

  /** set M to pointer newPtr
   *  \param newM the new M matrix
   */
  inline void setMPtr(SP::OSNSMatrix newM)
  {
    _M = newM;
  }

  // --- Q ---
  /** get the value of q, the constant vector in the LinearOSNS
   *  \return SiconosVector
   */
  inline const SiconosVector getQ() const
  {
    return *_q;
  }

  /** get q, the the constant vector in the LinearOSNS
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q;
  }

  /** set q to pointer newPtr
   *  \param newQ the new q vector
   */
  inline void setQPtr(SP::SiconosVector newQ)
  {
    _q = newQ;
  }

  /** get the type of storage for M
   * \return the type of storage
   */
  inline int getMStorageType() const
  {
    return _MStorageType;
  };

  /** set which type of storage will be used for M
   * \warning this function does not allocate any memory for M,
   * it just sets an indicator for future use
   * \param i an integer
   */
  inline void setMStorageType(int i)
  {
    _MStorageType = i;
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
  void computeq(double time);

  /** pre-treatment for LinearOSNS
      \param time the current time
      \return true if succeeded
   */
  virtual bool preCompute(double time);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *  \param time the current time
   *  \return information about the solver convergence.
   */
  virtual int compute(double time) = 0;

  /** post-treatment for LinearOSNS
   */
  virtual void postCompute();

  /** print the data to the screen
   */
  virtual void display() const;

  /** set if if _w and _z vectors are initialized with
      previous values of Y and Lambda when a change occurs in problem
      size
    * \param val true if we keep the previous values
    */
  void setKeepLambdaAndYState(bool val)
  {
    _keepLambdaAndYState = val ;
  }

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();


};

#endif // LinearOSNS_H
