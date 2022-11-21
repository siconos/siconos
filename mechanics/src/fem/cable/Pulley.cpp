#include "Pulley.h"
#include "Ropeway.h"



Pulley::Pulley(const Pile &a_pile) : Support(a_pile)
{
	m_radiusP = 0.;
	m_TR = 0.;
    dy = 0;
}


Pulley::~Pulley()
{}

void Pulley::prepare(const Rope &a_rope) {}

void Pulley::prepare(const Pile &a_start, const Pile &a_end, double T)
{
  Point delta;
  delta.diff(a_start, a_end);
  m_radiusP = delta.norm() / 2.0;  
  m_p.add(a_start, a_end);
  m_p.mult(0.5);  
  dy = a_start.y - m_p.y;
  m_TR = T;
}

int Pulley::compute(int nb, vector<Point> &a_q, int q_offset) const
{ 	
	double pi2 = M_PI_2;
	if (dy <= 0) {
	  pi2 = -M_PI_2;
	}
	size_t n = a_q.size();
	for (int i = 0; i < nb && (i + q_offset)<n; i++) {
	  double theta = pi2 + ((2 * M_PI_2) / (nb - 1)) * i;

	  Point &p = a_q[i + q_offset];
      p.x = m_p.x + m_radiusP * std::cos(theta);
      p.y = m_p.y + m_radiusP * std::sin(theta);
	  p.z = m_p.z;
	}	
	return q_offset + nb-1;
}

void Pulley::compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T)
{
  double dx = a_p.x - m_p.x;
  double dy = a_p.y - m_p.y;
  double go = sqrt(dx * dx + dy * dy) - get_radius();
  double nx = dx / (go + get_radius());
  double ny = dy / (go + get_radius());
  if (go <= a_tol) {
    g = go;
    G.x = nx;
    G.y = ny;
    T.x = -ny;
    T.y = nx;
  }
}

bool Pulley::isContact(const Point &a_p, const double &a_tol) 
{ 
  double dx = a_p.x - m_p.x;
  double dy = a_p.y - m_p.y;
  double go = sqrt(dx * dx + dy * dy) - get_radius();  
  return (go > a_tol);    
}

const double & Pulley::get_radius() const
{
	return m_radiusP;
}

const Point & Pulley::get_center()
{
	return m_p;
}


double Pulley::get_L(const Ropeway & a_rope) const
{
	return m_radiusP * M_PI / (1 + m_TR / a_rope.get_meca0().get_EA());
}

void to_json(ojson & j, const Pulley & p)
{
	j["radius"] = p.m_radiusP;
	j["center"] = p.m_p;
	j["tension"] = p.m_TR;
}
