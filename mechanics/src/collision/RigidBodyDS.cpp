#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "RigidBodyDS.hpp"
#include "SiconosContactor.hpp"
RigidBodyDS::RigidBodyDS(SP::SiconosVector position,
                         SP::SiconosVector velocity,
                         double mass,
                         SP::SimpleMatrix inertia)
  : NewtonEulerDS(position, velocity, mass, inertia)
  , _contactors(std::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
}

RigidBodyDS::~RigidBodyDS()
{
}
