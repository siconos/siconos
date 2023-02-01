#include "Carriers.h"

Carriers::Carriers()
{
  m_rho = 0.0;
  m_d = 0.0;
  m_mass = 0.0;
}

Carriers::~Carriers() {}

const double &Carriers::get_rho() const { return m_rho; }

const double &Carriers::get_d_inter_vehicules() const { return m_d; }

void Carriers::from_json(const json &j)
{
  j.at("rho").get_to(m_rho);
  j.at("mass").get_to(m_mass);
  j.at("d").get_to(m_d);
}
