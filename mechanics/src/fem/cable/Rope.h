#pragma once
#include "Cable.h"
#include "Pile.h"

class Rope {
public:
  /** Build a rope span (a piece of cable between to pylons)
      \param a_pile1 first pylon supporting the rope
      \param a_pile2 second pylon
      \param a_tol tolerance used to perform (stop) Newton-Raphson iterations
      \param n_max max. number of iterations used in Newton-Raphson
  */
  Rope(const Pile &a_pile0, const Pile &a_pile1, double a_tol, int n_max);

  virtual ~Rope();

  /** Compute initial profile of the Rope using Catenary equations

      Update positions, internal forces and tension in the rope.

      \param[in] a_meca mechanical properties of the rope span
      \param[in] nb_nodes number of nodes uses to compute profile with catenary equations
      \param[in] a_T0 initial tension applied at the beginning of the rope
      \para[in] m a_R0 support reaction at the beginning of the rope
  */
  void compute(const class Cable &a_meca, int nb_nodes, double a_T0, const Point &a_R0);

  int computeNbNodes(int nb_elem, double L);
  int computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R, std::vector<double> &a_TS,
                  int q_offset, bool a_reverse = false);

  double get_T0();
  double get_LastT();
  Point get_LastR();
  double get_L();
  const Point &get_SR() const;
  const Cable &get_meca() const;
  const Pile &get_pile0() const;

  void to_json(ojson &j);

private:
  Point ropeway_inc;
  std::vector<Point> q;   // positions
  std::vector<Point> R;   // internal forces [x,y,z]-> [H,V,B]
  std::vector<double> TS; // tension

  Cable meca;        // contient T0, EA, rho
  const Pile &pile0; // référence vers le support associé
  const Pile &pile1; // référence vers le support associé
  Point SR;          // support reaction [H,V,B]
  bool m_last;

  /** level of tolerance for the cable equation */
  double tol = 1e-7;
  /** max number of iterations */
  int n_max = 50;
  /** */
  int m_nbNodes;

  /**
      Modifies cable_inc such that the cable equation is satisfied (using tolerance tol and for
     a max number of iteration equal n_max) n_max iteration

      \param a_meca object which handles the geometric and material description of the cable
      \param bc a vector of two Pile objects
      \return a Point
  */
  Point get_adm_1C(const Cable &a_meca, const std::vector<Pile> &bc);

  /**
     \returns a Point, initial guess for the length and the slopes
     \param a vector of two Pile objects
   */
  Point guess(const std::vector<Pile> &bc);

  /**
     computes the residual equation and the jacobian of a fixed-fixed cable with imposed
     initial tension

     \param a_meca object which handles the geometric and material description of the cable
     \param bc a vector of two Pile objects
     \param cable_inc = [L, etaY, etaZ], with L
     the unstretched length of the cable, etaY the initial y-slope and etaZ the initial z-slope
     \param J the resulting jacobian
   */
  void cable_eq(const Cable &a_meca, const std::vector<Pile> &bc, const Point &cable_inc,
                Point &r, std::vector<std::vector<double>> &J);

  /** Computes initial profile of the cable

     \param[in] a__meca object which handles the geometric and material description of the
     cable \param[in] cable_inc [L, etaY, etaZ], with L the unstretched length of the cable,
     etaY the initial y-slope and etaZ the initial z-slope \param[in] nb_nodes number of nodes
     in the system \param[in, out] a_q positions \param[in, out] a_R internal forces
     \param[in,out] a_TS tension
     \param[in] q_offset
     \param[in] a_reverse

  */
  void get_profile_1C(const Cable &a_meca, const Point &cable_inc, int nb_nodes,
                      std::vector<Point> &a_q, std::vector<Point> &a_R,
                      std::vector<double> &a_TS, int q_offset = 0, bool a_reverse = false);
};
