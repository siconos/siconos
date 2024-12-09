#include "RigidBodyDS.hpp"

#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "SiconosContactor.hpp"
RigidBodyDS::RigidBodyDS(SP::SiconosVector position, SP::SiconosVector velocity, double mass,
                         SP::SimpleMatrix inertia)
    : NewtonEulerDS(position, velocity, mass, inertia),
      _contactors(std::make_shared<SiconosContactorSet>()),
      _useContactorInertia(true),
      _allowSelfCollide(true),
      _qExtrapolated(nullptr) {}

RigidBodyDS::~RigidBodyDS() {}

void RigidBodyDS::compute_extrapolated_position(double extrapolationCoefficient) {
  // we compute an extrapolation of the position
  if (!_qExtrapolated) _qExtrapolated.reset(new SiconosVector(q()->size()));

  SP::SiconosVector velocityIncrement(new SiconosVector(_twist->size()));

  _qExtrapolated->setValue(0, velocityIncrement->getValue(0));
  _qExtrapolated->setValue(1, velocityIncrement->getValue(1));
  _qExtrapolated->setValue(2, velocityIncrement->getValue(2));
  quaternionFromTwistVector(*velocityIncrement, *_qExtrapolated);
  const SiconosVector& qold = qMemory().getSiconosVector(0);
  compositionLawLieGroup(qold, *_qExtrapolated);
}
