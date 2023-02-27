#include "Support.h"

#include "TransportCable.h" // for sgn function

#ifndef NSICONOS
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothLaw.hpp"
#include "SiconosVector.hpp"
#endif

Support::Support(const Pile &a_pile) : r_pile(a_pile)
{
  m_p = r_pile;

#ifndef NSICONOS
  m_nslaw = nullptr;

  m_pc2 = std::make_shared<SiconosVector>(3);
  m_normal = std::make_shared<SiconosVector>(3);
  m_tangent = std::make_shared<SiconosVector>(3);
#endif
}

Support::~Support() {}

const double &Support::get_radius() const { return r_pile.get_radius(); }

void Support::prepare(const Rope &a_rope)
{
  double radius = sgn(a_rope.get_SR().z) * get_radius();
  m_p.z -= radius;
}

void Support::prepare(const Pile &a_start, const Pile &a_end, double T)
{
  // par défaut, ne fait rien
}

void Support::compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c)
{
  c = isContact(a_tol, a_p.x - m_p.x, 0, a_p.z - m_p.z, g, G.x, G.y, G.z, T.x, T.y, T.z) ? 1
                                                                                         : 0;
}

void to_json(ojson &j, const Support &s)
{
  j["radius"] = s.get_radius();
  j["p"] = s.m_p;
}

#ifndef NSICONOS
bool Support::isContact(const std::shared_ptr<SiconosVector> &a_p, const double &a_tol)
{
  double dx = a_p->getValue(0) - m_p.x;
  double dz = a_p->getValue(2) - m_p.z;
  double d = sqrt(dx * dx + dz * dz);
  double go = d - get_radius();
  bool isCt = (go <= a_tol);
  if (isCt) {
    m_normal->setValue(0, dx / d);
    m_normal->setValue(2, dz / d);
    m_tangent->setValue(0, -m_normal->getValue(2));
    m_tangent->setValue(2, m_normal->getValue(0));
  }
  return isCt;
}

void Support::InitFriction(double a_mu)
{
  if (m_nslaw != nullptr) {
    // m_nslaw->
  }
  else {
    m_nslaw = std::make_shared<NewtonImpactFrictionNSL>(0.0, 0.0, a_mu, 2);
  }
}

std::shared_ptr<NonSmoothLaw> Support::nslaw() { return m_nslaw; }

#endif

bool Support::isContact(const double &a_tol, const double &dx, const double &dy,
                        const double &dz, double &g, double &nx, double &ny, double &nz,
                        double &tx, double &ty, double &tz)
{
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  double go = d - get_radius();
  bool isCt = (go <= a_tol);
  if (isCt) {
    g = go;
    nx = dx / d;
    ny = dy / d;
    nz = dz / d;
    tx = -ny - nz;
    ty = (ny) ? nx : 0.0;
    tz = (nz) ? nx : 0.0;
  }
  return isCt;
}
