
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

void SiconosContactor::addShape(SP::SiconosPlane plane,
                                SP::SiconosVector offset)
{
  if (!offset) { offset.reset(new SiconosVector(7)); offset->zero(); }
  _planes.push_back(std::make_pair(plane, offset));
}

void SiconosContactor::addShape(SP::SiconosSphere sphere,
                                SP::SiconosVector offset)
{
  if (!offset) { offset.reset(new SiconosVector(7)); offset->zero(); }
  _spheres.push_back(std::make_pair(sphere, offset));
}

void SiconosContactor::addShape(SP::SiconosBox box,
                                SP::SiconosVector offset)
{
  if (!offset) { offset.reset(new SiconosVector(7)); offset->zero(); }
  _boxes.push_back(std::make_pair(box, offset));
}

void SiconosContactor::addShape(SP::SiconosConvexHull ch,
                                SP::SiconosVector offset)
{
  if (!offset) { offset.reset(new SiconosVector(7)); offset->zero(); }
  _chs.push_back(std::make_pair(ch, offset));
}
