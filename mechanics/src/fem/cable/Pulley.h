#pragma once
#include "Point.h"

class Pulley
{
public:
	Pulley();
	virtual ~Pulley();
	void init(const double &a_radius, const Point &a_center);

	const double &get_radius();
	const Point &get_center();
	

	void set_T(double T);
	double get_L(const class Ropeway &a_rope);

	friend void to_json(ojson &j, const Pulley &p);

private:
	Point center;
	double m_radiusP;
	double m_TR; // tension
};

