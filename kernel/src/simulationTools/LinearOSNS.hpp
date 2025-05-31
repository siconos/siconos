/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "SiconosVector.hpp"
#include "OneStepNSProblem.hpp"

// namespace siconos::numerics {

#include "NM_types.h"  // For NM_DENSE, NM_types ...
// enum NM_types : int;  // explicit type specification is required
//}  // namespace siconos::numerics


namespace siconos::modeling {
class NonSmoothLaw;
}

namespace siconos::nonsmooth_formulations {

class OSNSMatrix;

enum class LinearOSNSAssemblyType { REDUCED_BLOCK, GLOBAL, REDUCED_DIRECT, GLOBAL_REDUCED };

/**
   Base (abstract) class for linear non-smooth problems

   Usually in the form:

   \f$ w =  q + M z \f$

   where
   - \f$ w \in R^{n} \f$  and \f$ z \in R^{n} \f$ are the unknowns,
   - \f$ M \in R^{n \times n } \f$  and \f$ q \in R^{n} \f$

   examples: LCP, FrictionContact ...

*/
class LinearOSNS : public OneStepNSProblem {
 protected:
  ACCEPT_SERIALIZATION(LinearOSNS);

  /** vector w of a LinearOSNS system */
  std::shared_ptr<siconos::algebra::SiconosVector> _w{nullptr};

  /** vector z of a LinearOSNS system */
  std::shared_ptr<siconos::algebra::SiconosVector> _z{nullptr};

  /** matrix M of a LinearOSNS system */
  std::shared_ptr<OSNSMatrix> _M{nullptr};

  /** vector q of a LinearOSNS system */
  std::shared_ptr<siconos::algebra::SiconosVector> _q{nullptr};

  /** matrix W of a LinearOSNS system */
  std::shared_ptr<OSNSMatrix> _W{nullptr};

  /** matrix W of a LinearOSNS system */
  std::shared_ptr<OSNSMatrix> _W_inverse{nullptr};

  /** matrix H of a LinearOSNS system */
  std::shared_ptr<OSNSMatrix> _H{nullptr};

  /** Assembly strategy */
  LinearOSNSAssemblyType _assemblyType{LinearOSNSAssemblyType::REDUCED_BLOCK};

  /** Storage type for M - NM_DENSE: SiconosMatrix (dense), NM_SPARSE_BLOCK:
     Sparse Storage (embedded into OSNSMatrix) */
  NM_types _numericsMatrixStorageType{NM_DENSE};

  /** a boolean to decide if _w and _z vectors are initialized with
      previous values of Y and Lambda when a change occurs in problem
      size */
  bool _keepLambdaAndYState = true;

  /** nslaw effects : visitors experimentation
   */
  struct _TimeSteppingNSLEffect;
  struct _EventDrivenNSLEffect;
  struct _NSLEffectOnSim;
  friend struct _TimeSteppingNSLEffect;
  friend struct _EventDrivenNSLEffect;
  friend struct _NSLEffectOnSim;

  LinearOSNS() = default;

 public:
  /** constructor from a pre-defined solver options set.
   *
   *  \param options the options set
   *  \param assemblytype the method used to build the assembled matrix (enum
   * LinearOSNSAssemblyType)
   */
  LinearOSNS(std::shared_ptr<SolverOptions> options,
             LinearOSNSAssemblyType assemblyType = LinearOSNSAssemblyType::REDUCED_BLOCK)
      : OneStepNSProblem{options}, _assemblyType{assemblyType} {};

  /** destructor */
  virtual ~LinearOSNS() noexcept = default;

  /**
      current w vector (pointer link)

      \return pointer on a siconos::algebra::SiconosVector
  */
  inline std::shared_ptr<siconos::algebra::SiconosVector> w() const { return _w; }

  /**
     set w vector (pointer link)

     \param newPtr the new std::shared_ptr<siconos::algebra::SiconosVector>
  */
  inline void setWPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _w = newPtr; }

  /**
     current z vector (pointer link)

     \return pointer on a siconos::algebra::SiconosVector
  */
  inline std::shared_ptr<siconos::algebra::SiconosVector> z() const { return _z; }

  /**
     set z vector (pointer link)

     \param newPtr the new std::shared_ptr<siconos::algebra::SiconosVector>
  */
  inline void setzPtr(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) { _z = newPtr; }

  /**
      M matrix (pointer link)

      \return pointer on a OSNSMatrix
  */
  inline std::shared_ptr<OSNSMatrix> M() const { return _M; }

  /**
      set M to pointer newPtr

      \param newM the new M matrix
  */
  inline void setMPtr(std::shared_ptr<OSNSMatrix> newM) { _M = newM; }

  // --- H ---

  /** get H
   *
   *  \return pointer on a OSNSMatrix
   */
  inline std::shared_ptr<OSNSMatrix> H() const { return _H; }

  /** set the value of H
   *
   *  \param H the new matrix
   */
  void setH(std::shared_ptr<OSNSMatrix> H) { _H = H; }

  /**
     get q, the the constant vector in the LinearOSNS

     \return pointer on a siconos::algebra::SiconosVector
  */
  inline std::shared_ptr<siconos::algebra::SiconosVector> q() const { return _q; }

  /**
     set q to pointer newPtr

     \param newQ the new q vector
  */
  inline void setQPtr(std::shared_ptr<siconos::algebra::SiconosVector> newQ) { _q = newQ; }

  /**
     get the type of storage used for M

     \return NM_types (NM_DENSE, NM_SPARSE_BLOCK)
  */
  inline NM_types getMStorageType() const { return _numericsMatrixStorageType; };

  /** set which type of storage will be used for M
   *  \warning this function does not allocate any memory for M,
   *  it just sets an indicator for future use
   *
   *  \param i (NM_DENSE, NM_SPARSE_BLOCK)
   */
  inline void setMStorageType(NM_types i) { _numericsMatrixStorageType = i; };

  /** set which type of assembly will be used for M
   */
  inline void setAssemblyType(LinearOSNSAssemblyType assemblyType) {
    _assemblyType = assemblyType;
  };

  /** Memory allocation or resizing for z,w,q */
  void initVectorsMemory();

  /** initialize the _M matrix */
  virtual void initOSNSMatrix();

  /**
     To initialize the LinearOSNS problem(computes topology ...)

     \param sim the simulation owning this OSNSPB
  */
  void initialize(std::shared_ptr<siconos::simulation::Simulation> sim) override;

  /** compute extra-diagonal interactionBlock-matrix
   *
   *  \param ed an edge descriptor
   */
  void computeInteractionBlock(
      const siconos::graphs::InteractionsGraph::EDescriptor &ed) override;

  /** compute diagonal Interaction block
   *
   *  \param vd a vertex descriptor
   */
  void computeDiagonalInteractionBlock(
      const siconos::graphs::InteractionsGraph::VDescriptor &vd) override;

  /** compute matrix M */
  virtual void computeM();

  /** To compute a part of the q vector of the OSNS
   *
   *  \param vertex vertex (interaction) which corresponds to the considered block
   *  \param pos the position of the first element of yOut to be set
   */
  virtual void computeqBlock(siconos::graphs::InteractionsGraph::VDescriptor &vertex,
                             unsigned int pos);

  /** compute vector q
   *
   *  \param time the current time
   */
  virtual void computeq(double time);

  /**
     build problem coefficients (if required)

     \param time the current time
     \return true if the indexSet is not empty
  */
  bool preCompute(double time) override;

  /** update interactions variables (y and lambda) according to current problem
   *  found solutions.
   */
  void postCompute() override;

  /** print the data to the screen
   */
  void display() const override;

  /**
      choose initialisation behavior for w and z.

      \param val true: init w and z with previous values
      of y and lambda saved in interactions, false: init to 0.
  */
  void setKeepLambdaAndYState(bool val) { _keepLambdaAndYState = val; }

  virtual bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw &nslaw) = 0;
};
}  // namespace siconos::nonsmooth_formulations
#endif  // LinearOSNS_H
