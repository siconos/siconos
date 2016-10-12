
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

void SiconosContactor::addShape(SP::SiconosPlane plane,
                                SP::SiconosVector offset)
{
  _planes.push_back(std::make_pair(plane, offset));
}

void SiconosContactor::addShape(SP::SiconosSphere sphere,
                                SP::SiconosVector offset)
{
  _spheres.push_back(std::make_pair(sphere, offset));
}

void SiconosContactor::addShape(SP::SiconosBox box,
                                SP::SiconosVector offset)
{
  _boxes.push_back(std::make_pair(box, offset));
}

void SiconosContactor::addShape(SP::SiconosConvexHull ch,
                                SP::SiconosVector offset)
{
  _chs.push_back(std::make_pair(ch, offset));
}
