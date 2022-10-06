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
			//P.computeDS();
                  auto cable =
                      std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass);



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
