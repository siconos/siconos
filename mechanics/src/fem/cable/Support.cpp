#include "Support.h"



Support::Support(const Rope &a_rope)
	: m_pile(a_rope.get_pile0())
{
	m_p = m_pile;
	double radius = sgn(a_rope.get_SR().z) * get_radius();
	m_p.z -= radius;
}


Support::~Support()
{
}

const double & Support::get_radius() const
{
	return m_pile.get_radius();
}

void Support::compute(const Point &a_p, int k, double a_tol, double mu_s,
	double &g,
	int &act,
	Point &G,
	Point &T,
	Point &blocked,
	Point &blocked_value,
	double &eye_c)
{
	double dx = a_p.x - m_p.x;
	double dz = a_p.z - m_p.z;
	double go = sqrt(dx*dx + dz * dz) - get_radius();
	double nx = dx / (go+get_radius());
	double nz = dz / (go+get_radius());
	if (go <= a_tol) {
		g = go;
		act = k;
		eye_c = mu_s;
		G.x = nx;
		G.z = nz;
		T.x = -nz;
		T.z = nx;
		blocked.y = 1;
		blocked_value.y = m_p.y;
	}
}

void to_json(ojson & j, const Support & s)
{
	j["radius"] = s.get_radius();
	j["p"] = s.m_p;	
}
