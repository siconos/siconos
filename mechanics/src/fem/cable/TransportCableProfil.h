#pragma once
#include "TransportCableModel.h"
#include "TransportCableResult.h"
#include "Ropeway.h"
#include "Pulleywrap.h"
#include "Support.h"

class TransportCableProfil
{
public:
	TransportCableProfil(const TransportCableModel &a_model, TransportCableResult &a_results);
	virtual ~TransportCableProfil();

	void computeInitialProfil(int nb_nodes, 
		int nodes_per_pulley=0, 
		double a_tol = 1e-20,
		int a_nmax = 20);

	void computeFEM(int nb_elem = 1400, 
		double a_eps = 0.1,
		double a_tol = 1e-3, 
		double mu_s = 0.8, 
		double mu_p = 1.1);
	
private:
	const TransportCableModel &r_model;
	TransportCableResult &r_results;
	
	Pulleywrap puller12;
	Pulleywrap puller21;

	void compute_punct_load(int nb_elem, double Lc, double d_prop=0.8);
	void compute_ineq_constraint(const vector<Point> &a_X, double a_tol=1e-3, double mu_s=0.8, double mu_p=1.1);
};

