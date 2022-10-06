#pragma once
#include "TransportCableModel.h"
#include "TransportCableResult.h"
#include "Pulleywrap.h"

template <typename T> T getParam(const json &a_arg, const string &a_name, T a_default)
{
	T vRet = a_default;
	if (!a_arg.is_null()) {
		if (a_arg.contains(a_name)) {
			a_arg.at(a_name).get_to(vRet);
		}
	}	
	return vRet;
}

class TransportCableManager
{
public:
	TransportCableManager();
	virtual ~TransportCableManager();

	int importModel(const json& a_input, const string &a_filename="");
	int computeFEM(const json &a_args, const string &a_outfile, ojson &output);
	int exportTC(const json & a_args, const string & a_outfile, ojson & output);

	int simulation(const json & a_model, const json & a_args, const string &a_filename, const string & a_outfile, ojson & output);

private:
	TransportCableModel m_model;
	TransportCableResult m_results;

	void computeDS();

    void compute_mass(double elem_length, double elem_rho);
    void compute_external_load(double elem_length, double elem_rho);
};

