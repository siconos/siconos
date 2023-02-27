#pragma once
#include "Pile.h"
#include "Rope.h"

#ifndef NSICONOS
class SiconosVector;
class NonSmoothLaw;
#endif

class Support {
public:
  Support(const Pile &a_pile);
  virtual ~Support();
  virtual const double &get_radius() const;

  //------------ statique -------------
  virtual void prepare(const Rope &a_rope);
  virtual void prepare(const Pile &a_start, const Pile &a_end, double T);

  virtual void compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c);

#ifndef NSICONOS
  //------------ dynamique -------------
  virtual bool isContact(const std::shared_ptr<SiconosVector> &a_p, const double &a_tol);

  virtual void InitFriction(double a_mu);
  std::shared_ptr<NonSmoothLaw> nslaw();

  inline std::shared_ptr<SiconosVector> pc2() const { return m_pc2; }
  inline std::shared_ptr<SiconosVector> normal() const { return m_normal; }
  inline std::shared_ptr<SiconosVector> tangent() const { return m_tangent; }
#endif

  //------ Export ----------
  friend void to_json(ojson &j, const Support &p);

  bool isContact(const double &a_tol, const double &dx, const double &dy, const double &dz,
                 double &g,

                 double &nx, double &ny, double &nz,

                 double &tx, double &ty, double &tz);

protected:
  const Pile &r_pile; // reference to the cable model
  Point m_p;

#ifndef NSICONOS
  std::shared_ptr<NonSmoothLaw> m_nslaw;

  std::shared_ptr<SiconosVector> m_pc2;
  std::shared_ptr<SiconosVector> m_normal;
  std::shared_ptr<SiconosVector> m_tangent;
#endif
};
