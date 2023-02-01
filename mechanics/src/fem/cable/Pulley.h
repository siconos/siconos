#pragma once
#include "Point.h"
#include "Support.h"

class Pulley : public Support {
public:
  Pulley(const Pile &a_pile);
  virtual ~Pulley();
  virtual const double &get_radius() const;
  const Point &get_center();
  double get_L(const class Ropeway &a_rope) const;

  //------------ statique -------------
  virtual void prepare(const Rope &a_rope);
  virtual void prepare(const Pile &a_start, const Pile &a_end, double T);

  int compute(int nb, std::vector<Point> &a_q, int q_offset = 0) const;
  virtual void compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c);

#ifndef NSICONOS
  //------------ dynamique -------------
  virtual bool isContact(const std::shared_ptr<SiconosVector> &a_p, const double &a_tol);

#endif

  //------ Export ----------
  friend void to_json(ojson &j, const Pulley &p);

private:
  double m_radiusP;
  double m_TR; // tension
  double dy;
};
