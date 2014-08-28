#include "OccBody.hpp"
#include "OccBody_impl.hpp"

#include <gp_Vec.hxx>
#include <gp_Quaternion.hxx>

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

void OccBody::updateContactShapes()
{
  for (ContactShapes::iterator csi = _contactShapes->begin();
       csi != _contactShapes->end(); ++ csi)
  {
    OccContactShape& cs = **csi;

    SiconosVector& q = *_q;

    const gp_Vec& translat = gp_Vec(q(0), q(1), q(2));

    const gp_Quaternion& rota = gp_Quaternion(q(4), q(5), q(6), q(3));

    gp_Trsf transfo;

    transfo.SetTranslation(translat);
    transfo.SetRotation(rota);

    cs.Move(transfo);

    // cf code from Olivier
    //const TopLoc_Location& aLoc = cs.Location();
    //const gp_Trsf& T = aLoc.Transformation();
    //TopLoc_Location aLocWithoutList(T);
    //cs.Location(aLocWithoutList(T));

  }
}


const ContactShapes& OccBody::contactShapes() const
{
  return *_contactShapes;
};
