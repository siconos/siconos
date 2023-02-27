#pragma once
#include "Support.h"

class Ropeway {
public:
  Ropeway();
  virtual ~Ropeway();

  /** Initialize/compute profile of rope segments and save them

     \param a_meca object which handles the geometric and material description of the cable
     \param a_piles the vector of pylons, must be ordered!
     \param nb_nodes number of nodes by segment
     \param a_tol tolerance used to compute initial profile
     \apram a_nmax max number of iterations used to compute initial profile
  */
  void compute(const Cable &a_meca, const std::vector<Pile> &a_piles, int nb_nodes,
               double a_tol = 1e-20, int a_nmax = 20);


  /**
     Create and prepare all supports, from the list of declared ropes
     (a support for each 'lower' pile of each rope).
     
     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void prepareSupport(std::vector<std::shared_ptr<Support>> &a_supports,
                      int &a_pulleyIdx) const;

  int computeNbNodes(int nb_elem, double L);
  int computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R, std::vector<double> &a_TS,
                  int q_offset);

  const Pile &get_FirstPile();
  const Pile &get_LastPile();
  double get_T0();
  double get_LastT();
  double get_L();
  const Cable &get_meca0() const;

  int to_json(ojson &j);
  void set_Down(bool a_value);

private:

  /** Create new support from the lower pile of a given rope
     
     \param[in] a_rope the considered Rope
     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void addSupport(const Rope &a_rope, std::vector<std::shared_ptr<Support>> &a_supports,
                  int &a_pulleyIdx) const;

  std::vector<Rope> m_ropes;
  bool m_down;
};
