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
/*! \file
  Time-Stepping simulation with projections on constraints
*/
#ifndef TIMESTEPPINGDIRECTPROJECTION_H
#define TIMESTEPPINGDIRECTPROJECTION_H

#include "TimeStepping.hpp"



/** \class TimeSteppingProjectOn constraints
 *  \brief Time-Stepping scheme with a direct projection onto the constraint 
 *  thanks to the GGL augmentation of the system
 *
 *  For details, have a look on
 *   Projected event-capturing time-stepping schemes for nonsmooth mechanical systems 
 *   with unilateral contact and Coulomb's friction
 *   Vincent Acary 
 *   Computer Methods in Applied Mechanics and Engineering, Elsevier,
 *   2013, 256, pp. 224-250
 *
 * 
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Aug 2010
 *  
 
 */
class TimeSteppingDirectProjection : public TimeStepping
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(TimeSteppingDirectProjection);


  /** level of IndexSet on which we project (default =1 (activated contact))
   */
  unsigned int _indexSetLevelForProjection;

  /** Number of iteration of projection
   */
  unsigned int _nbProjectionIteration;

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

  /** disabled or enabled projection (Debug Projection) */
  unsigned int _doProj;
  unsigned int _doOnlyProj;


public:

  /** Constructor with the time-discretisation.
  *  \param a pointer to a timeDiscretisation (linked to the model
  *  that owns this simulation)
     \param a one step integrator
     \param a one step non smooth problem for the velocity formulation
     \param a one step non smooth problem for the position formulation
  */
  TimeSteppingDirectProjection(SP::TimeDiscretisation td,
                                   SP::OneStepIntegrator osi,
                                   SP::OneStepNSProblem osnspb_velo,
                                   SP::OneStepNSProblem osnspb_pos,
                                   unsigned int _level = 1);

  virtual void initOSNS();

  /** default constructor
   */
  TimeSteppingDirectProjection() {};

  virtual ~TimeSteppingDirectProjection();




  virtual void updateWorldFromDS()
  {
    ;
  }

  /** get the Number of iteration of projection
   */
  inline unsigned int nbProjectionIteration()
  {
    return _nbProjectionIteration;
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

  inline void setDoProj(unsigned int v)
  {
    _doProj = v;
  }
  inline void setDoOnlyProj(unsigned int v)
  {
    _doOnlyProj = v;
  }

  /**
   */
  void advanceToEvent();


  void nextStep();

  /*
   */
  void computeCriteria(bool * runningProjection);

  void newtonSolve(double criterion, unsigned int maxStep);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

DEFINE_SPTR(TimeSteppingDirectProjection)

#endif // TIMESTEPPINGDIRECTPROJECTION_H





