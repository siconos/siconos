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
/*! \file
  LsodarOSI solver (from odepack)
*/

#ifndef LsodarOSI_H
#define LsodarOSI_H

#include <vector>

#include "OneStepIntegrator.hpp"
#include "SiconosConst.hpp"  // MACHINE_PREC

/**
   LsodarOSI solver (odepack)

   Many parameters are required as input/output for LSODAR. See the
   documentation of this function in externals/odepack/opkdmain.f to have a full
   description of these parameters.
   Most of them are read-only parameters (ie can not be set by user).

   Except:

   - jt: Jacobian type indicator (1 means a user-supplied full Jacobian, 2
   means an internally generated full Jacobian). Default = 2.
   - itol, rtol and atol

   ITOL   = an indicator for the type of error control.

   RTOL   = a relative error tolerance parameter, either a scalar or array of
   length NEQ.

   ATOL   = an absolute error tolerance parameter, either a scalar or an array
   of length NEQ. Input only.
*/
class LsodarOSI : public OneStepIntegrator {
 private:
  ACCEPT_SERIALIZATION(LsodarOSI);

  static constexpr double ATOL_DEFAULT = 100 * siconos::internal::MACHINE_PREC;
  static constexpr double RTOL_DEFAULT = 10 * siconos::internal::MACHINE_PREC;

  /** neq, ng, itol, itask, istate, iopt, lrw, liw, jt
   *  See opkdmain.f and lsodar routine for details on those variables.
   */
  std::vector<int> _intData = {0, 0, 1, 1, 0, 0, 2};

  /** relative tolerance */
  std::vector<double> rtol = {RTOL_DEFAULT};
  /** absolute tolerance */
  std::vector<double> atol = {ATOL_DEFAULT};

  /** real work array */
  std::vector<double> rwork = {};
  /** integer work array */
  std::vector<int> iwork = {};
  /** integer array used for output of root information */
  std::vector<int> jroot = {};
  /** temporary vector to save x values */
  SP::BlockVector _xWork{nullptr};

  SP::SiconosVector _xtmp{nullptr};
  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;

 public:
  enum LsodarOSI_ds_workVector_id { FREE, WORK_LENGTH };

  enum LsodarOSI_interaction_workVector_id { OSNSP_RHS, WORK_INTERACTION_LENGTH };

  enum LsodarOSI_interaction_workBlockVector_id { xfree, BLOCK_WORK_LENGTH };

  /** Lsodar counter : Number of steps taken for the problem so far. */
  static int count_NST;
  /** Number of RHS evaluations for the problem so far. */
  static int count_NFE;

  /** Default and only constructor */
  LsodarOSI();

  /** destructor
   */
  ~LsodarOSI() noexcept = default;

  /** \return int parameters for lsodar */
  inline const std::vector<int> intData() const { return _intData; }

  /** \return int parameter number i
   *
   *  \param i index number (starting from 0)
   */
  inline int intData(unsigned int i) const { return _intData[i]; }

  /** set _intData[i]
   *
   *  \param i index number (starting from 0)
   *  \param newValue the new value
   */
  inline void setIntData(unsigned int i, int newValue) { _intData[i] = newValue; }

  /** \return relative tolerance parameter for lsodar */
  inline const std::vector<double> &getRtol() const { return rtol; }

  /** \return absolute tolerance parameter for lsodar*/
  inline const std::vector<double> &getAtol() const { return atol; }

  /** \return the maximum number of steps for one call */
  inline int getMaxNstep() const { return iwork[5]; }

  /** \return real work vector parameter for lsodar */
  inline const std::vector<double> &getRwork() const { return rwork; }

  /** \return iwork */
  inline const std::vector<int> &getIwork() const { return iwork; }

  /** \return root information */
  inline const std::vector<int> &getJroot() const { return jroot; }

  /** set Jt value, Jacobian type indicator. Excerpts from the lsodar
   *  documentation. 1 means a user-supplied full (neq by neq) jacobian. 2 means
   *  an internally generated (difference quotient) full jacobian (using neq
   *  extra calls to f per df/dy value). 4 means a user-supplied banded
   * jacobian. 5 means an internally generated banded jacobian (using ml+mu+1
   * extra calls to f per df/dy evaluation). if jt = 1 or 4, the user must
   * supply a subroutine jac (the name is arbitrary) as described above under
   * jac. if jt = 2 or 5, a dummy argument can be used.
   *
   *  \param newJT new value for the jt parameter.
   */
  inline void setJT(int newJT) { _intData[6] = newJT; };

  /** set itol, rtol and atol (tolerance parameters for lsodar)
   *
   *  \param newItol itol value
   *  \param newRtol rtol value
   *  \param newAtol atol value
   */
  void setTol(int newItol, std::vector<double> &&newRtol, std::vector<double> &&newAtol);

  /** set itol, rtol and atol (scalar tolerance parameters for lsodar)
   *
   *  \param newItol itol value
   *  \param newRtol rtol value
   *  \param newAtol atol value
   */
  void setTol(int newItol, double newRtol, double newAtol);

  /** set the maximum number of steps for one call of Lsodar
   *
   *  \param maxNumberSteps the maximum number of steps
   */
  void setMaxNstep(int maxNumberSteps);

  /** set the minimum and maximum step sizes
   *
   *  \param minStep minimum step size
   *  \param maxStep maximum step size
   */
  void setMinMaxStepSizes(double minStep, double maxStep);

  /** set maximum method order
   *
   *  \param maxorderNonStiff maximum order for nonstiff methods
   *  \param maxorderStiff maximum order for stiff methods
   */
  void setMaxOrder(int maxorderNonStiff, int maxorderStiff);

  /** update doubleData and iwork memory size, when changes occur in _intData.
   */
  void updateData();

  /** fill xWork with a double
   *
   *  \param size size of x array
   *  \param array x array of double
   */
  void fillXWork(int *size, double *array);

  /** compute rhs(t) for all dynamical systems in the set
   *
   *  \param t current time of simulation
   */
  void computeRhs(double t);

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
   *
   *  \param t current time of simulation
   *  \param DSG0 the graph of DynamicalSystem
   */
  void computeJacobianRhs(double t, DynamicalSystemsGraph &DSG0);

  void f(int *sizeOfX, double *time, double *x, double *xdot);

  void g(int *nEq, double *time, double *x, int *ng, double *gOut);

  void jacobianfx(int *, double *, double *, int *, int *, double *, int *);

  /** initialization of the integrator
   */
  void initialize() override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   *
   *  \param t time of initialization
   *  \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   *
   *  \param inter the interaction
   *  \param interProp the properties on the graph
   *  \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction &inter,
                                           InteractionProperties &interProp,
                                           DynamicalSystemsGraph &DSG) override;

  /** get the number of index sets required for the simulation
   *
   *  \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 3; };

  /** integrate the system, between tinit and tend (->iout=true), with possible
   *  stop at tout (->iout=false)
   *
   *  \param tinit initial time
   *  \param tend end time
   *  \param tout real end time
   *  \param ioparam in-out parameter, input: 1 for first call, else 2. Output:
   *  2 if no root was found, else 3.
   */
  void integrate(double &tinit, double &tend, double &tout, int &ioparam) override;

  /** update the state of the DynamicalSystems attached to this Integrator
   *
   *  \param level level of interest for the dynamics
   */
  void updateState(const unsigned int level) override;

  void prepareNewtonIteration(double time) override { assert(0); };

  /** integrates the Interaction linked to this integrator, without taking
   *  non-smooth effects into account
   *
   *  \param vertex_descr descriptor vertex of the interaction graph
   *  \param osnsp pointer to OneStepNSProblem
   */
  void computeFreeOutput(InteractionsGraph::VDescriptor &vertex_descr,
                         OneStepNSProblem *osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth
   * problem
   */
  SiconosVector &osnsp_rhs(InteractionsGraph::VDescriptor &vertex_inter,
                           InteractionsGraph &indexSet) override {
    return *(*indexSet.properties(vertex_inter).workVectors)[LsodarOSI::OSNSP_RHS];
  };

  /** print the data to the screen
   */
  void display() override;

  /** Return current number of rhs call (for all lsodar-like OSIs!)
   *
   *  \return int
   */
  static int count_rhs_call() { return count_NFE; }

  /** Return the number of lsodar steps already done (for all lsodar-like OSIs!)
   *
   *  \return int
   */
  static int count_steps() { return count_NST; }

  ACCEPT_STD_VISITORS();
};

#endif  // LsodarOSI_H
