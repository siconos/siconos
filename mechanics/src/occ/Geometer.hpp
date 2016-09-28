#ifndef Geometer_hpp
#define Geometer_hpp

#include <SiconosVisitor.hpp>
#include <iostream>

struct Geometer : public SiconosVisitor
{
  SPC::OccContactShape base;
  SP::ContactShapeDistance answer;

  Geometer() {};

  Geometer(const OccContactShape& base) : base(createSPtrConstOccContactShape(base)) {};

  using SiconosVisitor::visit;

  void visit(const OccContactFace& face)
  {
    answer = base->distance(face);
  }

  void visit(const OccContactEdge& edge)
  {
    answer = base->distance(edge);
  }

  virtual SP::ContactShapeDistance distance(const OccContactShape& psh1,
                                            const OccContactShape& psh2)
  {
    return psh1.distance(psh2);
  }
};


#endif
