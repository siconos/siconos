#ifndef MBTB_JOINTR
#define MBTB_JOINTR
#include "SiconosKernel.hpp"
#include "MechanicsFwd.hpp"
/**
 * \brief This class implements a joint in a multi-bodies system.
 * It is an aggregation to the class siconos::NewtonEulerR. Mainly, it consists in adding members needed for the computation of joint forces.
 */
class MBTB_JointR
{
  friend class PivotJointR;

protected:

public:
  MBTB_JointR();
  //! it is assumed that _interaction
  SP::Interaction  _interaction;
  //! it is assumed that _joinrR is a pivot.
  SP::NewtonEulerJointR _jointR;
  //! The first dynamical systems of the joint.
  SP::NewtonEulerDS _ds1;

  //!vector G0Ci, where Ci Contact points where the forces of joint must be computed.
  SP::SiconosVector _G0C1;
  //!vector G0Ci, where Ci Contact points where the forces of joint must be computed.
  SP::SiconosVector _G0C2;
  //! Joint forces. F1: F[0,1,2]. F2: F[3,4,5].
  SP::SiconosVector _F;
  //! A matrix such that  _M * F = BLambda
  SP::SimpleMatrix _M;


  //!It consists in building the system  _M * F = BLambda  and solving it.
  void computeEquivalentForces();




};


#endif
