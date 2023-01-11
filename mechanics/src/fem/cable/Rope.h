#pragma once
#include "Cable.h"
#include "Pile.h"

class Rope
{
public:
	Rope(const Pile &a_pile0, const Pile &a_pile1, double a_tol, int n_max);
	virtual ~Rope();

	void compute(const class Cable &a_meca,
		int nb_nodes,
		double a_T0,
		const Point &a_R0);
	
	int computeNbNodes(int nb_elem, double L);
	int computeMesh(vector<Point> &a_q, vector<Point> &a_R, vector<double> &a_TS, int q_offset, bool a_reverse = false);
	
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
	std::vector<Point> q;	// positions
	std::vector<Point> R;	// internal forces [x,y,z]-> [H,V,B]
	std::vector<double> TS; // tension

	Cable meca;			// contient T0, EA, rho
	const Pile &pile0;	// référence vers le support associé
	const Pile &pile1;	// référence vers le support associé
	Point SR;			// support reaction [H,V,B]
	bool m_last;
	
	double tol;
	int n_max;
	int m_nbNodes;

	Point get_adm_1C(const Cable &a_meca,
		const std::vector<Pile> &bc);

	Point guess(const std::vector<Pile> &bc);

	void cable_eq(const Cable &a_meca,
		const std::vector<Pile> &bc,
		const Point &cable_inc,
		Point &r,
		std::vector<std::vector<double>> &J);

	void get_profile_1C(const Cable &a_meca,
		const Point &cable_inc,
		int nb_nodes,
		vector<Point> &a_q,
		vector<Point> &a_R, 
		vector<double> &a_TS,
		int q_offset=0,
		bool a_reverse=false);
};

