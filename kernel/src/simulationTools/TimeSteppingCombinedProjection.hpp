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
/*! \file
  Time-Stepping simulation with projections on constraints
*/
#ifndef TIMESTEPPINGCOMBINEDPROJECTION_H
#define TIMESTEPPINGCOMBINEDPROJECTION_H

#include "TimeStepping.hpp"



/** Time-Stepping scheme
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.4.0.
 *  \date (Creation) May 2012
 *
 */
class TimeSteppingCombinedProjection : public TimeStepping
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(TimeSteppingCombinedProjection);

  /** level of IndexSet on which we project
   *  (default =2 (subset of activated constraint with positive reactions))
   */
  unsigned int _indexSetLevelForProjection;

  /** Cumulated Number of steps perfomed is the Newton Loop */
  unsigned int _cumulatedNewtonNbIterations;

  /** Number of iteration of projection
  */
  unsigned int _nbProjectionIteration;

  /** Number of cumulated iteration of projection
  */
  unsigned int _nbCumulatedProjectionIteration;

  /** Number of iteration for stabilizating indexsets
  */
  unsigned int _nbIndexSetsIteration;



  /** tolerance for the violation of the equality
   *  constraints at the  position level.
   */
  double _constraintTol;

  /** tolerance for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _constraintTolUnilateral;


  /** maximum violation for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _maxViolationUnilateral;

  /** maximum violation for the violation of the equality
   *  constraints at the  position level.
   */
  double _maxViolationEquality;




  /** Default maximum number of projection iteration*/
  unsigned int _projectionMaxIteration;

  /** Default maximum number of index set activation iteration*/
  unsigned int _kIndexSetMax;

  /** disabled or enabled projection (Debug Projection) */
  bool _doCombinedProj;

  /** disabled or enabled projection On Equality (or Unilateral) for unilateral constraints */
  bool _doCombinedProjOnEquality;

  /** Boolean to check if the index sets are stabilized in the Combined Projection Algorithm */
  bool _isIndexSetsStable;

  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param level unsigned int: the level of the set to be updated
   */
  void updateIndexSet(unsigned int level);


  struct _SimulationEffectOnOSNSP;
  friend struct _SimulationEffectOnOSNSP;


public:

  virtual void initOSNS();




  /** Constructor with the time-discretisation.
   *  \param td a pointer to a timeDiscretisation (linked to the model
   *  that owns this simulation)
   *  \param osi a one step integrator
   * \param osnspb_velo a one step non smooth problem for the velocity formulation
   *  \param osnspb_pos a one step non smooth problem for the position formulation
   *  \param _level
   */
  TimeSteppingCombinedProjection(SP::TimeDiscretisation td,
                                 SP::OneStepIntegrator osi,
                                 SP::OneStepNSProblem osnspb_velo,
                                 SP::OneStepNSProblem osnspb_pos,
                                 unsigned int _level = 2);


  /** default constructor
   */
  TimeSteppingCombinedProjection() {};

  virtual ~TimeSteppingCombinedProjection();

  virtual void updateWorldFromDS()
  {
    ;
  }
  /** get the Number of iteration of projection
   * \return unsigned int nbProjectionIteration
   */
  inline unsigned int nbProjectionIteration()
  {
    return _nbProjectionIteration;
  }
  /** get the Number of cumulated iteration of projection
   * \return unsigned int
   */
  inline unsigned int nbCumulatedProjectionIteration()
  {
    return _nbCumulatedProjectionIteration;
  }

  /** get the  Cumulated Number of steps perfomed in the Newton Loop
   * \return unsigned int
   */
  inline unsigned int cumulatedNewtonNbIterations()
  {
    return _cumulatedNewtonNbIterations;
  }

  /** get the Number of iteration for stabilizating indexsets
   * \return unsigned int
   */
  inline unsigned int nbIndexSetsIteration()
  {
    return _nbIndexSetsIteration;
  }

  inline void setConstraintTol(double v)
  {
    _constraintTol = v;
  }

  inline void setConstraintTolUnilateral(double v)
  {
    _constraintTolUnilateral = v;
  }

  inline double maxViolationUnilateral()
  {
    return _maxViolationUnilateral;
  }
  inline double maxViolationEquality()
  {
    return _maxViolationEquality;
  }


  inline void setProjectionMaxIteration(unsigned int v)
  {
    _projectionMaxIteration = v;
  }

  inline void setDoCombinedProj(unsigned int v)
  {
    _doCombinedProj = v;
  }

  inline bool doCombinedProjOnEquality()
  {
    return _doCombinedProjOnEquality;
  }


  /**
   */
  void advanceToEvent();
  /**
   */
  void advanceToEventOLD();

  /*
   */
  void computeCriteria(bool * runningProjection);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

DEFINE_SPTR(TimeSteppingCombinedProjection)

#endif // TIMESTEPPINGCOMBINEDPROJECTION_H





