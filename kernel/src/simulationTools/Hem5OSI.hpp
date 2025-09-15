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
/** \file Hem5OSI.hpp
 * Hem5OSI solver (from E. Hairer software lists)
 */
#ifndef Hem5OSI_H
#define Hem5OSI_H

#include <vector>

#include "OneStepIntegrator.hpp"
#include "SiconosConst.hpp"  // MACHINE_PREC
#include "SiconosFortran.h"  // for siconos::hairer

/** Hem5OSI solver (odepack)
 *
 *
 * Many parameters are required as input/output for LSODAR. See the
 * documentation of this function in externals/odepack/opkdmain.f to have a full
 * description of these parameters.  \n Most of them are read-only parameters
 * (ie can not be set by user). \n Except: \n
 *  - jt: Jacobian type indicator (1 means a user-supplied full Jacobian, 2
 * means an internally generated full Jacobian). \n Default = 2.
 *  - itol, rtol and atol \n
 *    ITOL   = an indicator for the type of error control. \n
 *    RTOL   = a relative error tolerance parameter, either a scalar or array of
 * length NEQ. \n ATOL   = an absolute error tolerance parameter, either a
 * scalar or an array of length NEQ.  Input only.
 */
class Hem5OSI : public OneStepIntegrator {
 private:
  ACCEPT_SERIALIZATION(Hem5OSI);
  static constexpr auto HEM5_ATOL_DEFAULT = 100 * siconos::internal::MACHINE_PREC;
  static constexpr auto HEM5_RTOL_DEFAULT = 10 * siconos::internal::MACHINE_PREC;
  static constexpr double INITIAL_GUESS_TS = 1.e-3;

  /** vector of integer data for the integrator
   * _intData[0] NQ size of the position vector q
   * _intData[1] NV size of the velocity  vector v NQ >= NQ
   * _intData[2] NU size of the external dynamic vector u
   * _intData[3] NL size of the Lagrange multiplier vector lambda
   * _intData[4] ITOL indicates whether RTOL and ATOL are scalar (ITOL=0), or
   * array of dimension NQ + NV + NU (ITOL=1) _intData[5] IOUT selects the dense
   * output formula _intData[6] LWK length of real array rwork _intData[7] LIWK
   * length of integer array iwork See hem5.f
   */
  std::vector<int> _intData = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  int _idid{0};

  /** relative tolerance */
  std::vector<double> rtol = {};
  /** absolute tolerance */
  std::vector<double> atol = {};
  /** real work array */
  std::vector<double> rwork = {};
  /** int work array */
  std::vector<int> iwork = {};

  double _timeStep{INITIAL_GUESS_TS};  // initial step size guess

  /** temporary vector to save q values */
  SP::BlockVector _qWork;
  /** temporary vector to save v values */
  SP::BlockVector _vWork;
  /** temporary vector to save v values */
  SP::BlockVector _uWork;
  /** temporary vector to save a values */
  SP::BlockVector _aWork;
  /** temporary vector to save lambda values */
  SP::BlockVector _lambdaWork;
  /** temporary vector to save forces values */
  SP::BlockVector _forcesWork;

  SP::SiconosVector _qtmp;
  SP::SiconosVector _vtmp;
  SP::SiconosVector _utmp;
  SP::SiconosVector _atmp;
  SP::SiconosVector _lambdatmp;
  SP::SiconosVector _forcestmp;

  /** nslaw effects
   */
  struct _NSLEffectOnFreeOutput;
  friend struct _NSLEffectOnFreeOutput;

  class Hem5OSI_impl {
   public:
    std::shared_ptr<Hem5OSI> hem5osi{nullptr};
    void fprob(int* IFCN, int* NQ, int* NV, int* NU, int* NL, int* LDG, int* LDF, int* LDA,
               int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR, int* INDGC, int* INDFLR,
               int* INDFLC, double* time, double* q, double* v, double* u, double* xl,
               double* G, double* GQ, double* F, double* GQQ, double* GT, double* FL,
               double* QDOT, double* UDOT, double* AM);
    void solout(int* MODE, int* NSTEP, int* NQ, int* NV, int* NU, int* NL, int* LDG, int* LDF,
                int* LDA, int* LRDO, int* LIDO, siconos::fortran::hairer::fprobpointer FPROB,
                double* q, double* v, double* u, double* DOWK, int* IDOWK);

    Hem5OSI_impl(std::shared_ptr<Hem5OSI> h) : hem5osi(h) {}
  };

 public:
  std::shared_ptr<Hem5OSI_impl> _impl{nullptr};

  enum Hem5OSI_ds_workVector_id { FREE, WORK_LENGTH };

  enum Hem5OSI_interaction_workVector_id { OSNSP_RHS, WORK_INTERACTION_LENGTH };

  enum Hem5OSI_interaction_workBlockVector_id { xfree, BLOCK_WORK_LENGTH };

  /** Default and only constructor */
  Hem5OSI();

  /** destructor
   */
  ~Hem5OSI() noexcept = default;

  /** get vector of int parameters for lsodar
   *  \return a vector<int>
   */
  inline const std::vector<int> intData() const { return _intData; }

  /** get _intData[i]
   * \param i index
   * \return an int
   */
  inline int intData(unsigned int i) const { return _intData[i]; }
  /** set _intData[i]
   * \param i index
   * \param newValue
   */
  inline void setIntData(unsigned int i, int newValue) { _intData[i] = newValue; }

  /** \return relative tolerance parameter for Hem5 */
  inline const std::vector<double>& getRtol() const { return rtol; }

  /** \return absolute tolerance parameter for Hem5 */
  inline const std::vector<double>& getAtol() const { return atol; }

  /** get the maximum number of steps for one call
   *\return an interger
   */
  inline int getMaxNstep() const { return iwork[11]; }

  /** \return real work vector parameter for lsodar */
  inline const std::vector<double>& getRwork() const { return rwork; }

  /** \return iwork */
  inline const std::vector<int>& getIwork() const { return iwork; }

  /** set itol, rtol and atol (tolerance parameters for Hem5)
   *  \param itol int (itol value)
   *  \param rtol double * (rtol)
   *  \param atol double * (atol)
   */
  void setTol(int itol, std::vector<double>&& rtol, std::vector<double>&& atol);

  /** set itol, rtol and atol (scalar tolerance parameters for Hem5)
   *  \param itol int (itol value)
   *  \param rtol double (rtol)
   *  \param atol double (atol)
   */
  void setTol(int itol, double rtol, double atol);

  /** set the maximul number of steps for one call of Hem5OSI
   *\param nstepmax an int
   */
  void setMaxNstep(int nstepmax);

  /** set the minimum and maximum step sizes
   * \param maxstepsize double (maximul step size)
   */
  void setMaxStepSize(double maxstepsize);

  /** update _intData
   */
  void updateIntData();

  /** update doubleData and iwork memory size, when changes occur in _intData.
   */
  void updateData();

  /** fill qWork with a double
   *  \param sizex int*, size of x array
   *  \param x double* x:array of double
   */
  void fillqWork(int* sizex, double* x);

  /** fill vWork with a double
   *  \param sizex int*, size of x array
   *  \param x double* x:array of double
   */
  void fillvWork(int* sizex, double* x);

  /** compute rhs(t) for all dynamical systems in the set
   */
  void computeRhs(double);

  /** compute jacobian of the rhs at time t for all dynamical systems in the set
   */
  void computeJacobianRhs(double);

  unsigned int numberOfConstraints();

  void f(int* sizeOfX, double* time, double* x, double* xdot);

  void g(int* nEq, double* time, double* x, int* ng, double* gOut);

  void jacobianfx(int*, double*, double*, int*, int*, double*, int*);

  /** initialization of the integrator
   */
  void initialize() override;
  /** initialization of the work vectors and matrices (properties) related to
   *  one dynamical system on the graph and needed by the osi
   * \param t time of initialization
   * \param ds the dynamical system
   */
  void initializeWorkVectorsForDS(double t, SP::DynamicalSystem ds) override;

  /** initialization of the work vectors and matrices (properties) related to
   *  one interaction on the graph and needed by the osi
   * \param inter the interaction
   * \param interProp the properties on the graph
   * \param DSG the dynamical systems graph
   */
  void initializeWorkVectorsForInteraction(Interaction& inter,
                                           InteractionProperties& interProp,
                                           DynamicalSystemsGraph& DSG) override;

  /** get the number of index sets required for the simulation
   * \return unsigned int
   */
  unsigned int numberOfIndexSets() const override { return 3; };

  /** integrate the system, between tinit and tend (->iout=true), with possible
   * stop at tout (->iout=false) \param tinit initial time \param tend end time
   *  \param tout real end time
   *  \param idid in-out parameter, input: 1 for first call, else 2. Output: 2
   * if no root was found, else 3.
   */
  void integrate(double& tinit, double& tend, double& tout, int& idid) override;

  /** update the state of the DynamicalSystems attached to this Integrator
   *  \param level level of interest for the dynamics
   */
  void updateState(const unsigned int level) override;

  void prepareNewtonIteration(double time) override { assert(0); };

  /** integrates the Interaction linked to this integrator, without taking
   *  non-smooth effects into account
   *
   *  \param vertex_inter of the interaction graph
   *  \param osnsp pointer to OneStepNSProblem
   */
  void computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter,
                         OneStepNSProblem* osnsp) override;

  /** return the workVector corresponding to the right hand side of the OneStepNonsmooth
   * problem
   */
  SiconosVector& osnsp_rhs(InteractionsGraph::VDescriptor& vertex_inter,
                           InteractionsGraph& indexSet) override {
    return *(*indexSet.properties(vertex_inter).workVectors)[Hem5OSI::OSNSP_RHS];
  };

  /** print the data to the screen
   */
  void display() override;

  ACCEPT_STD_VISITORS();
};

#endif  // Hem5OSI_H
