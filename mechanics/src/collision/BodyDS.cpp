
#include "MechanicsFwd.hpp"
#include "NewtonEulerDS.hpp"
#include "BodyDS.hpp"
#include "SiconosContactor.hpp"
#include "SiconosContactorBase.hpp"

#include <boost/make_shared.hpp>

BodyDS::BodyDS(SP::SiconosVector position,
               SP::SiconosVector velocity,
               double mass,
               SP::SimpleMatrix inertia)
  : NewtonEulerDS(position, velocity, mass, inertia)
  , SiconosContactorBase(SP::SiconosVector())
  , _contactors(std11::make_shared<SiconosContactorSet>())
{
}

BodyDS::~BodyDS()
{
}
