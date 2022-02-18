/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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


/*! \file DiskMovingPlanR.hpp
 */
/** \class DiskMovingPlanR
 *  \brief disk - moving plan relation - Inherits from LagrangianRheonomousR
 */

#ifndef DiskMovingPlanR_h
#define DiskMovingPlanR_h

#include "MechanicsFwd.hpp"
#include "LagrangianRheonomousR.hpp"
#include "PluggedObject.hpp"

typedef double(*FTime)(double);


#define COMPUTE(X) \
  { if (_##X##Function->fPtr) _##X=((FTime)(_##X##Function->fPtr))(t); else _##X=0.; }



class DiskMovingPlanR : public LagrangianRheonomousR,
  public std::enable_shared_from_this<DiskMovingPlanR>
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(DiskMovingPlanR);

  double _time, _A, _B, _C, _ADot, _BDot, _CDot, _sqrA2pB2, _r, _AADot, _BBDot, _cubsqrA2pB2;


  SP::PluggedObject _AFunction;
  SP::PluggedObject _BFunction;
  SP::PluggedObject _CFunction;

  SP::PluggedObject _ADotFunction;
  SP::PluggedObject _BDotFunction;
  SP::PluggedObject _CDotFunction;

  DiskMovingPlanR() : LagrangianRheonomousR() {};

public:

  DiskMovingPlanR(FTime, FTime, FTime, FTime, FTime, FTime, double);

  void init(double);

  using LagrangianRheonomousR::computeh;
  /** to compute the output y = h(t,q,z) of the Relation
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
      \param y the resulting vector
  */
  void computeh(double time, const BlockVector& q, BlockVector& z, SiconosVector& y) override;

  /** to compute the jacobian of h(...). Set attribute _jachq (access: jacqhq())
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  void computeJachq(double time, const BlockVector& q, BlockVector& z) override;

  using LagrangianRheonomousR::computehDot;
  /** to compute the time-derivative of the output y = h(t,q,z), saved in attribute _hDot (access: hDot())
      \param time current time value
      \param q coordinates of the dynamical systems involved in the relation
      \param z user defined parameters (optional)
  */
  void computehDot(double time, const BlockVector& q, BlockVector& z) override;

  double distance(double, double, double);

  void setComputeAFunction(FTime f)
  {
    _AFunction.reset(new PluggedObject());
    _AFunction->setComputeFunction((void*) f);
  }

  void setComputeBFunction(FTime f)
  {
    _BFunction.reset(new PluggedObject());
    _BFunction->setComputeFunction((void*) f);
  }

  void setComputeCFunction(FTime f)
  {
    _CFunction.reset(new PluggedObject());
    _CFunction->setComputeFunction((void*) f);
  }

  void setComputeADotFunction(FTime f)
  {
    _ADotFunction.reset(new PluggedObject());
    _ADotFunction->setComputeFunction((void*) f);
  }

  void setComputeBDotFunction(FTime f)
  {
    _BDotFunction.reset(new PluggedObject());
    _BDotFunction->setComputeFunction((void*) f);
  }

  void setComputeCDotFunction(FTime f)
  {
    _CDotFunction.reset(new PluggedObject());
    _CDotFunction->setComputeFunction((void*) f);
  }

  bool equal(FTime, FTime, FTime, double) const;

  /** compute A
      \param t the time
  */
  void computeA(double t)
  COMPUTE(A)

  /** compute B
      \param t the time
  */
  void computeB(double t)
  COMPUTE(B)

  /** compute C
      \param t the time
  */
  void computeC(double t)
  COMPUTE(C)
    
  /** compute ADot
    \param t the time
  */
  inline void computeADot(double t)
  COMPUTE(ADot)

  /** compute BDot
      \param t the time
  */
  inline void computeBDot(double t)
  COMPUTE(BDot)


  /** compute CDot
      \param t the time
  */
  inline void computeCDot(double t)
  COMPUTE(CDot)

  /** visitor hooks
   */
  ACCEPT_VISITORS();

  ~DiskMovingPlanR() {};

};
#undef COMPUTE

#endif /* DiskMovingPlanR */

