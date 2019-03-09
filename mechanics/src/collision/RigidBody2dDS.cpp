
#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "RigidBody2dDS.hpp"
#include "SiconosContactor.hpp"

#include <boost/make_shared.hpp>

RigidBody2dDS::RigidBody2dDS(SP::SiconosVector position,
                             SP::SiconosVector velocity,
                             SP::SiconosMatrix  mass)
  : LagrangianDS(position, velocity, mass)
  , _contactors(std11::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
  // Check size of positions, velocities and mass matrix
}

RigidBody2dDS::~RigidBody2dDS()
{
}
