#include "TransportCableProfil.h"

TransportCableProfil::TransportCableProfil(const TransportCableModel &a_model,
                                           TransportCableResult &a_results)
    : r_model(a_model), r_results(a_results)
{
}

void TransportCableProfil::computeInitialProfil(int nb_nodes, double a_tol, int a_nmax)
{
  // calcul les positions, tensions des cordes
  Cable meca = r_model.get_cable();
  Carriers vehicules = r_model.get_carriers();
  meca.set_rho(meca.get_rho() + vehicules.get_rho());

  r_results.rope1.compute(meca, r_model.get_piles1(), nb_nodes, a_tol, a_nmax);
  meca.set_T(r_results.rope1.get_T0());
  r_results.rope2.compute(meca, r_model.get_piles2(), nb_nodes, a_tol, a_nmax);

  // prépare les supports: pile, station -> support
  r_results.prepareSupport();
}

void TransportCableProfil::computeFEM(int nb_elem, double a_eps, double a_tol)
{
  std::shared_ptr<Pulley> puller12 =
      std::dynamic_pointer_cast<Pulley>(r_results.supports[r_results.puller12idx]);

  std::shared_ptr<Pulley> puller21 =
      std::dynamic_pointer_cast<Pulley>(r_results.supports[r_results.puller21idx]);

  double Lt = puller12->get_L(r_results.rope2);
  double Lb = puller21->get_L(r_results.rope1);
  r_results.length = r_results.rope1.get_L() + Lt + r_results.rope2.get_L() + Lb;

  int n_Pt = (int)rint(nb_elem * Lt / r_results.length);
  int n_Pb = (int)rint(nb_elem * Lb / r_results.length);

  int n1 = r_results.rope1.computeNbNodes(nb_elem, r_results.length);
  int n2 = r_results.rope2.computeNbNodes(nb_elem, r_results.length);
  r_results.nb_nodes = n_Pt + n1 + n_Pb + n2;

  std::vector<Point> &q = r_results.q;
  q.clear();
  q.resize(r_results.nb_nodes);
  std::vector<Point> &R = r_results.R;
  R.clear();
  R.resize(r_results.nb_nodes);
  std::vector<double> &TS = r_results.TS;
  TS.clear();
  TS.resize(r_results.nb_nodes);
  int offset = r_results.rope1.computeMesh(q, R, TS, 0);
  offset = puller12->compute(n_Pt + 1, q, offset);
  offset = r_results.rope2.computeMesh(q, R, TS, offset);
  puller21->compute(n_Pb + 1, q, offset);

  r_results.elem_length = r_results.length / r_results.nb_nodes;
  compute_punct_load(r_results.nb_nodes, r_results.length);

  compute_ineq_constraint(q, a_tol);
  /*dgi = np.zeros(nb_node)
  dgi[gi <= 0] = gi[gi <= 0]
  q = q - np.matmul(np.transpose(Gi), 1.1*dgi)*/
  std::vector<double> &g = r_results.g;
  std::vector<std::vector<Point>> &G = r_results.G;
  double k = 1 + a_eps;
  for (auto i = 0; i < r_results.nb_nodes; i++) {
    Point p;
    if (g[i] < 0) {
      for (auto j = 0; j < r_results.nb_nodes; j++) {
        p.add(p, G[j][i]);
      }
      p.mult(k * g[i]);
      q[i].diff(q[i], p);
    }
  }
}

void TransportCableProfil::compute_punct_load(int nb_elem, double Lc)
{
  /*
  get_punct_load(nb_elem,Lc,rho_vehicules,d_inter_vehicules)

  where

  nb_elem is the number of element in the mesh
  Lc is the total lentgh of cable
  rho_vehicules is the fictitious linear density due to vehicules and d_inter_vehicules
  is the distance between each vehicules (unstretched length)

  Returns a vector which n-th element is the mass hanged to node n the first vehicule is placed
  randomly between 0 and d_inter_vehicules on the cable
  */

    std::vector<double>& punct = r_results.punct;
    punct.clear();
    punct.resize(nb_elem, 0);

  const Carriers &vehicules = r_model.get_carriers();
  double d = vehicules.get_d_inter_vehicules();
  if (d) {
      std::vector<double> lc = linspace(0.0, Lc, nb_elem);

      int nb_vehicules = (int)(Lc / d);
      double m = vehicules.get_rho() * d;

      double d_prop = vehicules.get_d_start();
      double start = d;
      if (d_prop < 0 || d_prop > 1)
          start *= (double)rand() / RAND_MAX;
      else
          start *= d_prop;
      
      int si = 0;
      while (si < nb_elem && lc[si] < start) si++;
      if (si < nb_elem) {
          int k = 1;
          punct[si] = m;
          for (auto i = si; i < nb_elem; i++) {
              int ind = (int)((lc[i] - start) / d);
              if (ind > k) {
                  punct[i] = m;
                  k++;
                  if (k > nb_vehicules)
                      break;
              }
          }
      }
  }  
}

void TransportCableProfil::compute_ineq_constraint(const std::vector<Point> &a_X, double a_tol)
{
  /*
  @author: charl

  get_ineq_constraint(q,supports)

  where

  X is the coordinates of cable particles
  supports is the list of obstacles which i-th element is containing [ positions, radius ]
  associated to piles (cylinder with circled base in xy plane) pulleys is the list of pulleys
  which i-th element is containing [ positions, radius ] associated to pulleys (cylinder with
  circled base in xz plane)

  Return the inequality constraint vector for the support contained into supports and the
  active set vector which i-th element is :
          - NaN if the constraint is inactive
          - Obstacle index in supports if the constraint is active
  */
  r_results.prepareIneqConstraint((int)a_X.size());

  std::vector<int> &c = r_results.contacts;
  std::vector<double> &g = r_results.g;
  std::vector<std::vector<Point>> &G = r_results.G;
  std::vector<std::vector<Point>> &T = r_results.T;
  std::vector<std::shared_ptr<Support>> &supports = r_results.supports;

  for (auto &s : supports) {
    size_t i = 0;
    for (auto &p : a_X) {
      s->compute(p, a_tol, g[i], G[i][i], T[i][i], c[i]);
      i++;
    }
  }
}

int TransportCableProfil::to_json(ojson& ro)
{
  // j = {}
   
  // r_results.rope1.to_json(out);
  // r_results.rope2.to_json(out);
  return 0;
}
