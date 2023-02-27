#pragma once
#include "Pulley.h"
#include "Ropeway.h"
#include "Support.h"

class TransportCableResult {
public:
  TransportCableResult();
  
  virtual ~TransportCableResult() noexcept = default;

  void prepareSupport();
  void prepareIneqConstraint(int nb_nodes);

  int exportTC(const std::string &a_fileName, ojson &a_output,
               const std::string &a_option = "all");

  int to_json(ojson &j, const std::string &a_option = "all");

  int puller12idx{-1};
  int puller21idx{-1};
  Ropeway rope1;
  Ropeway rope2;

  std::vector<std::shared_ptr<Support>> supports;

  std::vector<Point> q;      // positions
  std::vector<Point> R;      // internal forces [x,y,z]-> [H,V,B]
  std::vector<double> TS;    // tension
  std::vector<int> contacts; // points en contact (=1)

  int nb_nodes{0};
  double length{0.};
  double elem_length{0.};
  // à convertir en siconos (vecteur ou matrice)
  std::vector<double> punct;

  std::vector<double> g;
  std::vector<std::vector<Point>> G;
  std::vector<std::vector<Point>> T;

#ifndef NSICONOS
  std::shared_ptr<class SiconosVector> q0{nullptr};
  std::shared_ptr<class SiconosVector> v0{nullptr};

  std::shared_ptr<class SiconosMatrix> mass{nullptr};
  std::shared_ptr<class SiconosVector> b{nullptr};
#endif
};
