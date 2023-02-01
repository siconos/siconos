#include "Cable.h"

Cable::Cable()
{
  m_rho = 0.0;
  m_EA = 1.0;
  m_T0 = 1.0;
}

Cable::~Cable() {}

const double &Cable::get_EA() const { return m_EA; }

const double &Cable::get_rho() const { return m_rho; }

const double &Cable::get_T0() const { return m_T0; }

double Cable::get_alpha() const { return 9.81 * m_rho / m_T0; }

double Cable::get_beta() const { return m_T0 / m_EA; }

void Cable::set_T(double a_T)
{
  if (a_T > 0.0)
    m_T0 = a_T;
}

void Cable::set_rho(double a_rho)
{
  if (a_rho > 0.0)
    m_rho = a_rho;
}

void Cable::from_json(const json &j)
{
  j.at("EA").get_to(m_EA);
  j.at("rho").get_to(m_rho);
  j.at("T0").get_to(m_T0);
}
