
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

void SiconosContactor::addShape(SP::SiconosShape shape,
                                SP::SiconosVector offset)
{
  _shapes.push_back(ShapeOffset(shape, offset));
}

// Needed for quaternion calculations below
#include <bullet/LinearMath/btQuaternion.h>
#include <bullet/LinearMath/btVector3.h>

void SiconosContactor::setPosition(const SP::SiconosVector position)
{
  std::vector<ShapeOffset>::iterator s;
  for (s = _shapes.begin();
       s != _shapes.end(); s++)
  {
    /* Adjust offset position according to current rotation */
    SP::SiconosVector pos(new SiconosVector(7));
    btQuaternion rbase((*position)(4), (*position)(5),
                       (*position)(6), (*position)(3));
    btVector3 rboffset = quatRotate(rbase, btVector3((*s->offset)(0),
                                                     (*s->offset)(1),
                                                     (*s->offset)(2)));

    /* Calculate total orientation */
    btQuaternion roffset((*s->offset)(4), (*s->offset)(5),
                         (*s->offset)(6), (*s->offset)(3));
    btQuaternion r(rbase * roffset);

    /* Set the absolute shape position */
    pos->setValue(0, position->getValue(0) + rboffset.x());
    pos->setValue(1, position->getValue(1) + rboffset.y());
    pos->setValue(2, position->getValue(2) + rboffset.z());
    pos->setValue(3, r.w());
    pos->setValue(4, r.x());
    pos->setValue(5, r.y());
    pos->setValue(6, r.z());

    s->shape->setPosition(pos);
  }
}
