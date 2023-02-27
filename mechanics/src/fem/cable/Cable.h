#pragma once
#include "BaseModel.h"

class Cable : public BaseModel {
public:
  Cable();
  virtual ~Cable();

  const double &get_EA() const;
  const double &get_rho() const;
  const double &get_T0() const;

  double get_alpha() const;
  double get_beta() const;

  void set_T(double a_T);
  void set_rho(double a_rho);

private:
  void from_json(const json &j);

  double m_EA;  // rigidity (N)
  double m_rho; // linear density	(kg/m)
  double m_T0;  // initial tension (kN)
};
