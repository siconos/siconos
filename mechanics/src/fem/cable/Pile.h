#pragma once
#include "Point.h"
class Pile : public Point {
public:
  Pile();
  Pile(const Pile &a_pile, bool a_isStation);
  virtual ~Pile();
  const double &get_radius() const;
  // const double &get_dUp() const;
  // const double &get_dDown() const;

  void from_json(const json &j);
  bool isStation() const;
  void transform(bool a_Up);

private:
  bool m_isStation;
  double m_radius; // rayon de la poulie dans le cas station, rayon courbure sinon
  double m_dUp;    // distance poteau - cable montant
  double m_dDown;  // distance poteau - cable descendant
  double m_h;      // hauteur du poteau
};
