#pragma once
#include "Rope.h"

class Support
{
public:
	Support(const Rope &a_rope);
	virtual ~Support();
	const double &get_radius() const;

	void compute(const Point &a_p, int k, double a_tol, double mu_s,
		double &g, 
		int &act, 
		Point &G,
		Point &T,
		Point &blocked,
		Point &blocked_value,
		double &eye_c);

	friend void to_json(ojson &j, const Support &p);

private:
	const Pile &m_pile;
	Point m_p;
};

