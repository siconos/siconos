#ifndef WHICHGEOMETER_HPP
#define WHICHGEOMETER_HPP

#include "MechanicsFwd.hpp"
#include "Question.hpp"
#include "Geometer.hpp"

template <typename T>
struct WhichGeometer : public Question<SP::Geometer>
{
  using SiconosVisitor::visit;

  void visit(const OccContactFace& s)
  {
    answer.reset(new FaceGeometer<T>(s));
  }

  void visit(const OccContactEdge& s)
  {
    answer.reset(new EdgeGeometer<T>(s));
  }
};

#endif
