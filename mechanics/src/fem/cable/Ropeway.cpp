#include "Ropeway.h"
#include "Cable.h"
#include <algorithm>


Ropeway::Ropeway()
{
	m_down = false;
}


Ropeway::~Ropeway()
{
}

void Ropeway::compute(const Cable &a_meca, const vector<Pile> &a_piles, int nb_nodes, double a_tol, int a_nmax)
{
	/* get_ropeway(meca, piles)

	where

	meca is a list containing [T0, EA, rho]
		T0  : initial tension
		EA  : rigidity of the cable
		rho : density of the cable

	piles is a np.array which i-th row contains a pile coordinate np.array([px, py, pz])
		This list should be ordered by the physical sequence of support

	nb_nodes is the number of nodes by subspan

	Returns the positions in the shape of a 3*nb_nodes vector [x,y,z,...,x,y,z],
	the tension in the shape of a nb_nodes vector [T,...T],
	the internal force vector in the shape of a 3*nb_nodes vector [H,V,B,...,H,V,B]
	the support reaction of piles in the shape of a 6*nb_piles vector [px,py,pz,H,V,B]
	Moreover, it provides with admissibilities for each subspan with a np.array which i-th row
	contains cable unknowns for each subspan np.array([L, etaY, etaZ])
*/
	size_t n = a_piles.size();
	m_ropes.clear();
	double T0 = a_meca.get_T0();
	Point R0;	
	for (size_t k = 0; k < n; k++) {		
		size_t k1 = k;
		if (k < n - 1) {
			k1++;
		}
		Rope r(a_piles[k], a_piles[k1], a_tol, a_nmax);
		r.compute(a_meca, nb_nodes, T0, R0);
		T0 = r.get_LastT();
		R0 = r.get_LastR();
		m_ropes.push_back(r);		
	}
}

void Ropeway::prepareSupport(vector<Support>& a_supports) const
{
	if(!m_down) {
		for (auto &r : m_ropes) {
			if (!r.get_pile0().isStation())
				a_supports.push_back(Support(r));
		}
	}
	else {
		for (auto r = m_ropes.rbegin(); r != m_ropes.rend(); r++) {
			if (!r->get_pile0().isStation())
				a_supports.push_back(Support(*r));
		}
	}
}

int Ropeway::computeNbNodes(int nb_elem, double L)
{
	int N = 0;
	for (auto &r : m_ropes) {
		N += r.computeNbNodes(nb_elem, L);
	}
	return N;
}

int Ropeway::computeMesh(vector<Point>& a_q, int q_offset)
{
	int offset = q_offset;
	if (!m_down) {
		for (auto &r : m_ropes) {
			offset += r.computeMesh(a_q, offset);
		}
	}
	else {
		for (auto r = m_ropes.rbegin(); r != m_ropes.rend(); r++) {			
			offset += r->computeMesh(a_q, offset, true);
		}
	}	
	return offset;
}

double Ropeway::get_T0()
{
	if (m_ropes.size()) {
		return m_ropes.front().get_T0();
	}
	else
		return 0.0;
}

double Ropeway::get_LastT()
{
	if (m_ropes.size()) {
		return m_ropes.back().get_LastT();
	}
	else
		return 0.0;
}

double Ropeway::get_L()
{
	double l = 0.0;
	for (auto &r : m_ropes) {
		l += r.get_L();
	}
	return l;
}

const Cable & Ropeway::get_meca0() const
{
	return m_ropes.front().get_meca();
}

int Ropeway::to_json(ojson & j)
{
	j = {
				{"ropeway_inc",  ojson::array()},
				{"q",  ojson::array()},
				{"TS",  ojson::array()},
				{"R",  ojson::array()},
				{"SR",  ojson::array()},
				{"meca_global",  ojson::array()}
	};
	for (auto &r : m_ropes) {
		r.to_json(j);
	}
	return 0;
}

void Ropeway::set_Down(bool a_value)
{
	m_down = a_value;
}
