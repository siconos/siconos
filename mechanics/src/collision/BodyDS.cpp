
#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "BodyDS.hpp"
#include "SiconosContactor.hpp"

#include <boost/make_shared.hpp>

BodyDS::BodyDS(SP::SiconosVector position,
               SP::SiconosVector velocity,
               double mass,
               SP::SimpleMatrix inertia)
  : NewtonEulerDS(position, velocity, mass, inertia)
  , _contactors(std11::make_shared<SiconosContactorSet>())
  , _useContactorInertia(true)
  , _allowSelfCollide(true)
{
}

BodyDS::~BodyDS()
{
}
