#pragma once
#include "TransportCable.h"
#include "Pulley.h"

class Pulleywrap
{
public:
	Pulleywrap(const Point &start, const Point &end, Pulley &a_pulley);
	virtual ~Pulleywrap();

	void compute(int nb);
	int compute(int nb, vector<Point> &a_q, int q_offset=0);

	void compute(const vector<Point> &a_p, int k, double a_tol, double mu,
		class TransportCableResult &a_result);


	const double &get_radius();
	void set_T(double T);
	double get_L(const class Ropeway &a_rope);

	int to_json(ojson &j);

private:	
	const Point &r_start;
	const Point &r_end;
	Pulley &r_pulley;

	std::vector<Point> qP;
		
	void compute(const Point &a_p, int k, double a_tol, double mu,
		double &g,
		int &act,
		Point &G,
		Point &T,
		Point &blocked,
		Point &blocked_value,
		double & eye_c);
};

