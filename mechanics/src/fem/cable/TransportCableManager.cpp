#include "TransportCableManager.h"
#include "TransportCableProfil.h"
#include "CableDS.hpp"
//#include "TCException.h"


TransportCableManager::TransportCableManager()
{
}

TransportCableManager::~TransportCableManager()
{
}

int TransportCableManager::importModel(const json & a_input, const string & a_filename)
{
	int res = EXIT_FAILURE;
	if (a_input.is_null()) {
		res = m_model.from_file(a_filename);
	}
	else {
		res = m_model.from_json(a_input);
	}
	return res;
}

int TransportCableManager::computeFEM(const json & a_args, const string & a_outfile, ojson & output)
{
	if (m_model.isLoaded()) {
        TransportCableProfil P(m_model, m_results);
        string method = "all";
		method = getParam(a_args, "compute", method);
		
		P.computeInitialProfil(
			getParam(a_args, "nb_node0", 50),
			getParam(a_args, "nodes_per_pulley", 0),
			getParam(a_args, "tol", 1e-20),
			getParam(a_args, "nmax", 20)
		);
		
		P.computeFEM(
			getParam(a_args, "nb_node", 1400),
			getParam(a_args, "eps", 0.1),
			getParam(a_args, "tol_contact", 1e-3),
			getParam(a_args, "mu_s", 0.8),
			getParam(a_args, "mu_p", 1.1)
		);

		if (method == "all") {
			computeDS();                 
		}
		exportTC(a_args, a_outfile, output);

		return EXIT_SUCCESS;
	}
	else {
		//throw TCException("Load a model before compute");		
	}
	return EXIT_FAILURE;
}

int TransportCableManager::exportTC(const json & a_args, const string & a_outfile, ojson & output)
{
	string vOption = getParam(a_args, "export", (string)"all");
	m_results.exportTC(a_outfile, output, vOption);

	return EXIT_SUCCESS;
}

int TransportCableManager::simulation(const json & a_model, const json & a_args, const string & a_filename, const string & a_outfile, ojson & output)
{
	int vRet = importModel(a_model, a_filename);
	if (vRet == EXIT_SUCCESS) {	
		vRet = computeFEM(a_args, a_outfile, output);
	}
	return vRet;
}

void TransportCableManager::computeDS()
{
	// model is loaded
	// q0 must be computed
	// q0 = q	
	int ndof = m_results.q.size()*3;
	m_results.q0 = std::make_shared<SiconosVector>(ndof); 
	m_results.v0 = std::make_shared<SiconosVector>(ndof); 
	size_t i = 0;
	for (auto &q : m_results.q) { 
		m_results.q0->setValue(i++, q.x);
        m_results.q0->setValue(i++, q.y);
        m_results.q0->setValue(i++, q.z);        
    }
		
	double rho = m_model.get_cable().get_rho();

	// compute mass -> results.mass
	compute_mass(m_results.elem_length, rho);

	// compute external forces -> results.b
    compute_external_load(m_results.elem_length, rho);

	// create dynamics model
    auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(m_results.q0, 
																	m_results.v0,
                                                                    m_results.mass, 
		m_model.get_cable().get_EA(),
        m_results.elem_length
		);

	cable->setFExtPtr(m_results.b);


}


void TransportCableManager::compute_mass(double a_length, double a_rho)
{
  /*get_mass_damp(elem_length,elem_rho,damping)
	where
	elem_length is the initial length of each element
	elem_rho is the vector which e-th element is the linear density of element number e

  */
	int ndof = m_results.q0->size();
	m_results.mass = std::make_shared<SimpleMatrix>(ndof, ndof, Siconos::UBLAS_TYPE::SPARSE);
	double k = a_rho * a_length / 3.0;
	for (size_t i = 0; i < ndof - 3; i++) {
		m_results.mass->setValue(i, i, 4 * k);
		m_results.mass->setValue(i, i + 3, k);
		m_results.mass->setValue(i + 3, i, k);
	}
	for (size_t i = 0; i < 3; i++) {
		m_results.mass->setValue(i + ndof - 3, i + ndof - 3, 4 * k);
		m_results.mass->setValue(i, i + ndof - 3, k);
		m_results.mass->setValue(i + ndof - 3, i, k);
	}
}

void TransportCableManager::compute_external_load(double a_length, double a_rho)
{
	size_t ndof = m_results.q0->size();
	m_results.b = std::make_shared<SiconosVector>(ndof);
  	  
	double k = -9.81 * a_rho * a_length;
	for (size_t i = 2; i < ndof; i += 3) {
         m_results.b->setValue(i, k);
	}
}

