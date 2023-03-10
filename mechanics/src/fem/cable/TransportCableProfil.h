#pragma once
#include "Support.h"
#include "TransportCableModel.h"
#include "TransportCableResult.h"

class TransportCableProfil {
public:
  TransportCableProfil(const TransportCableModel &a_model, TransportCableResult &a_results);
  
  virtual ~TransportCableProfil() noexcept = default;

  void computeInitialProfil(int nb_nodes, double a_tol = 1e-20, int a_nmax = 20);

  void computeFEM(int nb_elem, double a_eps = 0.1, double a_tol = 1e-3);

  int to_json(ojson& out);
  
private:
  const TransportCableModel &r_model;
  TransportCableResult &r_results;

  void compute_punct_load(int nb_elem, double Lc);

  /**
     \param a_X vector of positions, coordinates of cable 'particles'
     \param a_tol tolerance used to activate contacts
   */
  void compute_ineq_constraint(const std::vector<Point> &a_X, double a_tol = 1e-3);
};
