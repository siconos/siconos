#include "Pulleywrap.h"
#include "Ropeway.h"
#include "TransportCableResult.h"


Pulleywrap::Pulleywrap(const Point &start, const Point &end, Pulley &a_pulley)
	: r_start(start), r_end(end), r_pulley(a_pulley)
{	
}


Pulleywrap::~Pulleywrap()
{
}

void Pulleywrap::compute(int nb)
{
	qP.resize(nb);
	compute(nb, qP);
}

int Pulleywrap::compute(int nb, vector<Point> &a_q, int q_offset)
{
	/*
@author: charl

get_pulley_wrap(start,end)

where

start is the entrance point of the cable in the pulley

end is the exit point of the pulley

and nb is the number of nodes

Returns the cable profile qP, the radius of the pulley radiusP and the center of it.
Be careful, the profile is following the order given in input and start and end should have
the same horizontal and vertical components
*/	
	Point delta;	
	delta.diff(r_start, r_end);
	double vRadiusP = delta.norm() / 2.0;
	Point center;
	center.add(r_start, r_end);
	center.mult(0.5);
	double dy = r_start.y - center.y;
	double pi2 = M_PI_2;
	if (dy <= 0) {
		pi2 = -M_PI_2;
	}
	size_t n = a_q.size();
	for (int i = 0; i < nb && (i + q_offset)<n; i++) {
		double theta = pi2 + ((2 * M_PI_2) / (nb - 1)) * i;

		Point &p = a_q[i + q_offset];
		p.x = center.x + vRadiusP * std::cos(theta);		
		p.y = center.y + vRadiusP * std::sin(theta);
		p.z = center.z;
	}	
	r_pulley.init(vRadiusP, center);

	return q_offset + nb-1;
}

void Pulleywrap::compute(const vector<Point>& a_p, int k, double a_tol, double mu,
	TransportCableResult &a_result)
{
	vector<double>& g = a_result.g;
	vector<int>& act = a_result.act;
	vector<vector<Point>>& G = a_result.G;
	vector<vector<Point>>& T = a_result.T;
	vector<Point>& blocked = a_result.blocked;
	vector<Point>& blocked_value = a_result.blocked_value;
	vector<double>& eye_c = a_result.eye_c;
	size_t i = 0;
	for (auto &p : a_p) {
		compute(p, k, a_tol, mu, g[i], act[i], G[i][i], T[i][i], blocked[i], blocked_value[i], eye_c[i]);
		i++;
	}
}

void Pulleywrap::compute(const Point &a_p, int k, double a_tol, double mu,
	double &g,
	int &act,
	Point &G,
	Point &T,
	Point &blocked,
	Point &blocked_value,
	double & eye_c)
{
	const Point &center = r_pulley.get_center();
	double dx = a_p.x - center.x;
	double dy = a_p.y - center.y;
	double go = sqrt(dx*dx + dy * dy) - get_radius();
	double nx = dx / (go + get_radius());
	double ny = dy / (go + get_radius());
	if (go <= a_tol) {
		g = go;
		act = k;
		eye_c = mu;
		G.x = nx;
		G.y = ny;
		T.x = -ny;
		T.y = nx;
		blocked.z = 1;
		blocked_value.z = center.z;
	}
}

const double & Pulleywrap::get_radius()
{
	return r_pulley.get_radius();
}

void Pulleywrap::set_T(double T)
{
	r_pulley.set_T(T);
}

double Pulleywrap::get_L(const Ropeway & a_rope)
{
	return r_pulley.get_L(a_rope);
}

int Pulleywrap::to_json(ojson & j)
{
	j = {
				{"qP",  ojson::array()},
				{"radius",  get_radius()},
				{"center",  r_pulley.get_center()}
	};
	for (auto &q : qP) {
		j["qP"].push_back(q.x);
		j["qP"].push_back(q.y);
		j["qP"].push_back(q.z);
	}
	return EXIT_SUCCESS;
}
