#include "OccBody.hpp"
#include "OccBody_impl.hpp"

OccBody::OccBody(SP::SiconosVector position,
                 SP::SiconosVector velocity,
                 double mass ,
                 SP::SiconosMatrix inertia) :
  NewtonEulerDS(position, velocity, mass, inertia),
  _contactShapes(new ContactShapes())
  {};


void OccBody::addContactShape(SP::OccContactShape shape)
{
  _contactShapes->push_back(shape);
};

const ContactShapes& OccBody::contactShapes() const
{
  return *_contactShapes;
};
