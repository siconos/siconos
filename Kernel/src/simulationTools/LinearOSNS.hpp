/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file LinearOSNS.h
  \brief Linear Complementarity Problem formulation and solving
*/

#ifndef LinearOSNS_H
#define LinearOSNS_H

#include "OneStepNSProblem.hpp"
#include "SimpleVector.hpp"
#include "SimpleMatrix.hpp"
#include "BlockCSRMatrix.hpp"
#include <sys/time.h>

#include "SiconosPointers.hpp"

/** stl vector of double */
typedef std::vector<double> MuStorage;
TYPEDEF_SPTR(MuStorage);

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
  friend class _TimeSteppingNSLEffect;
  friend class _EventDrivenNSLEffect;
  friend class _NSLEffectOnSim;

  /** default constructor (private)
   */
  LinearOSNS() ;
  virtual void updateOSNSMatrix();
public:

  /** xml constructor
      \param SP::OneStepNSProblemXML : the XML linked-object
      \param string, ns problem type            \
  */
  LinearOSNS(SP::OneStepNSProblemXML, const std::string&);

  /** constructor from data
        \param int, numerics_solver_id
        \param string, ns problem type    \
    */
  LinearOSNS(const int ,
             const std::string&);
  /** constructor from data
      \param int, numerics_solver identifier
      \param string, ns problem type    \
      \param String: id of the problem (default = "unamed")
  */
  LinearOSNS(const int,
             const std::string&,
             const std::string& = "unamed");

  /** destructor
   */
  ~LinearOSNS() {};

  // --- W ---
  /** get the value of w, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an
   *  lvalue => return SimpleVector
   */
  inline const SimpleVector getW() const
  {
    return *_w;
  }

  /** get w, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SP::SiconosVector W() const
  {
    return _w;
  }

  /** set the value of w to newValue
   *  \param SiconosVector newValue
   */
  void setW(const SiconosVector&);

  /** set w to pointer newPtr
   *  \param SP::SiconosVector  newPtr
   */
  inline void setWPtr(SP::SiconosVector newPtr)
  {
    _w = newPtr;
  }

  // --- Z ---
  /** get the value of z, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an
   *  lvalue => return SimpleVector
   */
  inline const SimpleVector getz() const
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
   *  \param SiconosVector newValue
   */
  void setz(const SiconosVector&);

  /** set z to pointer newPtr
   *  \param SP::SiconosVector  newPtr
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

  /** set the value of M to newValue
   *  \param newValue OSNSMatrix
   */
  void setM(const OSNSMatrix&);

  /** set M to pointer newPtr
   *  \param newPtr OSNSMatrix*
   */
  inline void setMPtr(SP::OSNSMatrix newPtr)
  {
    _M = newPtr;
  }

  // --- Q ---
  /** get the value of q, the constant vector in the LinearOSNS
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an
   *  lvalue => return SimpleVector
   */
  inline const SimpleVector getQ() const
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
   *  \param SP::SiconosVector  newPtr
   */
  inline void setq(SP::SiconosVector newPtr)
  {
    _q = newPtr;
  }

  /** get the type of storage for M */
  inline int getMStorageType() const
  {
    return _MStorageType;
  };

  /** set which type of storage will be used for M - Note that this
      function does not allocate any memory for M, it just sets an
      indicator for future use */
  inline void setMStorageType(int i)
  {
    _MStorageType = i;
  };

  /** Memory allocation or resizing for z,w,q */
  void initVectorsMemory();

  /** Build or reinit M */
  virtual void updateM();

  /** To initialize the LinearOSNS problem(computes topology ...)
      \param the simulation, owner of this OSNSPB
  */
  virtual void initialize(SP::Simulation);

  /** compute extra-diagonal unitaryBlock-matrix
   *  \param an edge descriptor
   */
  virtual void computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor&);

  /** compute diagonal unitary block
   * \param a vertex descriptor
   */
  virtual void computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor&);

  /** To compute a part of the "q" vector of the OSNS
      \param SP::UnitaryRelation, the UR which corresponds to the considered block
      \param unsigned int, the position of the first element of yOut to be set
  */
  virtual void computeqBlock(SP::UnitaryRelation, unsigned int);

  /** compute vector q
   *  \param double : current time
   */
  void computeq(double);

  /** pre-treatment for LinearOSNS
   *  \param double : current time
   *  \return void
   */
  void preCompute(double);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
   *  \param double : current time
   *  \return int, information about the solver convergence.
   */
  virtual int compute(double) = 0;

  /** post-treatment for LinearOSNS
   */
  void postCompute() ;

  /** print the data to the screen
   */
  virtual void display() const;

  /** copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** set if if _w and _z vectors are initialized with
      previous values of Y and Lambda when a change occurs in problem
      size */
  void setKeepLambdaAndYState(bool val)
  {
    _keepLambdaAndYState = val ;
  }

};

TYPEDEF_SPTR(LinearOSNS);

#endif // LinearOSNS_H
