#pragma once
#include "BaseModel.h"

// punctual masses
class Carriers : public BaseModel {
public:
  Carriers();
  virtual ~Carriers();
  const double &get_rho() const;
  const double &get_d_inter_vehicules() const;
  const double& get_d_start() const;

private:
  void from_json(const json &j);

  int m_n;       // number of vehicules
  double m_mass; // mass of one vehicule (kg)
  double m_d;    // distance between two vehicules (m)
  double m_rho;

  double m_loaded_mass; // mass 100% of one vehicule (kg)
  double m_up_load;  // % up load
  double m_down_load; // % down load      
  double m_d_start; // % (m_d), distance of the first vehicule
};
