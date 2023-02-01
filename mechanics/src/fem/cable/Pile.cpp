#include "Pile.h"

Pile::Pile()
{
  m_radius = 0.;
  m_h = 0;
  m_dDown = 0;
  m_dUp = 0;
  m_isStation = false;
}

Pile::Pile(const Pile &a_pile, bool a_isStation) : Point(a_pile)
{
  m_radius = a_pile.get_radius();
  m_h = 0;
  m_dDown = 0;
  m_dUp = 0;
  m_isStation = a_isStation;
}

Pile::~Pile() {}

const double &Pile::get_radius() const { return m_radius; }

void Pile::from_json(const json &j)
{
  Point::from_json(j);
  j.at("R").get_to(m_radius);
  try {
    j.at("dUp").get_to(m_dUp);
    j.at("dDown").get_to(m_dDown);
    j.at("h").get_to(m_h);
  }
  catch (const std::exception &) {
  }
}

bool Pile::isStation() const { return m_isStation; }

void Pile::transform(bool a_Up)
{
  if (!a_Up) {
    if (m_isStation) {
      y += 2.0 * m_radius;
    }
    else {
      y += m_dUp + m_dDown;
    }
  }
}
