/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

  virtual void initOSNS();

  /** level of IndexSet on which we project
   *  (default =2 (subset of activated constraint with positive reactions))
   */
  unsigned int _indexSetLevelForProjection;

  /** tolerance for the violation of the equality
  *  constraints at the  position level.
  */
  double _constraintTol;

  /** tolerance for the violation of the unilateral
   *  constraints at the  position level.
   */
  double _constraintTolUnilateral;

  /** Default maximum number of projection iteration*/
  unsigned int _projectionMaxIteration;

  /** disabled or enabled projection (Debug Projection) */
  unsigned int _doCombinedProj;


  /** Boolean to check if the index sets are stabilized in the Combined Projection Algorithm */
  bool _isIndexSetsStable;

  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param unsigned int: the number of the set to be updated
   */
  void updateIndexSet(unsigned int);



public:

  /** Constructor with the time-discretisation.
  *  \param a pointer to a timeDiscretisation (linked to the model
  *  that owns this simulation)
     \param a one step integrator
     \param a one step non smooth problem for the velocity formulation
     \param a one step non smooth problem for the position formulation
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

  /** Update the _level* attributes and add the specific index Set
   * \param inter a new SP::Interaction
   * \param init bool to determine if we are in the initialisation phase
   */
  void computeLevelsForInputAndOutput(SP::Interaction inter, bool init = false);

  virtual void updateWorldFromDS()
  {
    ;
  }

  inline void setConstraintTol(double v)
  {
    _constraintTol = v;
  }

  inline void setConstraintTolUnilateral(double v)
  {
    _constraintTolUnilateral = v;
  }

  inline void setProjectionMaxIteration(unsigned int v)
  {
    _projectionMaxIteration = v;
  }

  inline void setDoCombinedProj(unsigned int v)
  {
    _doCombinedProj = v;
  }

  /**
   */
  void advanceToEvent();

  /*
   */
  void computeCriteria(bool * runningProjection);

  /** visitors hook
   */
  //ACCEPT_STD_VISITORS();

};

DEFINE_SPTR(TimeSteppingCombinedProjection);

#endif // TIMESTEPPINGCOMBINEDPROJECTION_H





