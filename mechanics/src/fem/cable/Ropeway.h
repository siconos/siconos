#pragma once
#include "Support.h"

class Ropeway
{
public:
	Ropeway();
	virtual ~Ropeway();

	void compute(
		const Cable &a_meca,
		const vector<Pile> &a_piles, 
		int nb_nodes, 
		double a_tol = 1e-20, 
		int a_nmax = 20);
        void prepareSupport(vector<std::shared_ptr<Support>> &a_supports,
                            int &a_pulleyIdx) const;

	int computeNbNodes(int nb_elem, double L);
	int computeMesh(vector<Point> &a_q, vector<Point> &a_R, vector<double> &a_TS, int q_offset);

	const Pile &get_FirstPile();
	const Pile &get_LastPile();
	double get_T0();
	double get_LastT();
	double get_L();
	const Cable &get_meca0() const;

	int to_json(ojson &j);
	void set_Down(bool a_value);

private:
        void addSupport(const Rope &a_rope, vector<std::shared_ptr<Support>> &a_supports,
                        int &a_pulleyIdx) const;
	vector<Rope> m_ropes;
	bool m_down;
};

