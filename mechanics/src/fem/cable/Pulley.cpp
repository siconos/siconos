#include "Pulley.h"
#include "Ropeway.h"



Pulley::Pulley()
{
	m_radiusP = 0.;
	m_TR = 0.;
}


Pulley::~Pulley()
{
}

void Pulley::init(const double & a_radius, const Point & a_center)
{
	m_radiusP = a_radius;
	center = a_center;
}

const double & Pulley::get_radius()
{
	return m_radiusP;
}

const Point & Pulley::get_center()
{
	return center;
}

void Pulley::set_T(double T)
{
	m_TR = T;
}

double Pulley::get_L(const Ropeway & a_rope)
{
	return m_radiusP * M_PI / (1 + m_TR / a_rope.get_meca0().get_EA());
}

void to_json(ojson & j, const Pulley & p)
{
	j["radius"] = p.m_radiusP;
	j["center"] = p.center;
	j["tension"] = p.m_TR;
}
