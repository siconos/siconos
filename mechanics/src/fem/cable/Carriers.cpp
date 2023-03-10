#include "Carriers.h"

Carriers::Carriers()
{
	m_n = 0;
  m_rho = 0.0;
  m_d = 0.0;
  m_mass = 0.0;

  m_loaded_mass = 0.0;
  m_up_load = 1;
  m_down_load = 0;
  m_d_start = 0;
}

Carriers::~Carriers() {}

const double &Carriers::get_rho() const { return m_rho; }

const double &Carriers::get_d_inter_vehicules() const { return m_d; }

const double& Carriers::get_d_start() const
{
	return m_d_start;
}

void Carriers::from_json(const json &j)
{
  j.at("rho").get_to(m_rho);
  j.at("mass").get_to(m_mass);
  j.at("d").get_to(m_d);

  if (j.contains("d_start")) {
	  j.at("d_start").get_to(m_d_start);
  }
}
