#include "Support.h"



Support::Support(const Pile &a_pile) : r_pile(a_pile)
{
	m_p = r_pile;	
}


Support::~Support()
{
}

const double & Support::get_radius() const
{
	return r_pile.get_radius(); }


void Support::prepare(const Rope &a_rope)
{
  double radius = sgn(a_rope.get_SR().z) * get_radius();
  m_p.z -= radius;
}

void Support::prepare(const Pile &a_start, const Pile &a_end, double T) {}


void Support::compute(
	const Point &a_p, 
	double a_tol, 
	double &g,
	Point &G,
	Point &T)
{
	double dx = a_p.x - m_p.x;
	double dz = a_p.z - m_p.z;
	double go = sqrt(dx*dx + dz * dz) - get_radius();
	double nx = dx / (go+get_radius());
	double nz = dz / (go+get_radius());
	if (go <= a_tol) {
		g = go;	
		G.x = nx;
		G.z = nz;
		T.x = -nz;
		T.z = nx;		
	}
}

bool Support::isContact(const Point &a_p, const double &a_tol) 
{ 
	double dx = a_p.x - m_p.x;
	double dz = a_p.z - m_p.z;
	double go = sqrt(dx * dx + dz * dz) - get_radius();  
	return (go > a_tol);		
}

void Support::getContact(const Point &a_p, SiconosVector &pc2, SiconosVector &normal,
                         SiconosVector &tangent)
{
  double dx = a_p.x - m_p.x;
  double dz = a_p.z - m_p.z;
  double go = sqrt(dx * dx + dz * dz) - get_radius();
  normal.setValue(0, dx / (go + get_radius()));
  normal.setValue(2, dz / (go + get_radius()));
  tangent.setValue(0, -normal.getValue(2));
  tangent.setValue(2, normal.getValue(0));
  
}

void to_json(ojson & j, const Support & s)
{
	j["radius"] = s.get_radius();
	j["p"] = s.m_p;	
}
