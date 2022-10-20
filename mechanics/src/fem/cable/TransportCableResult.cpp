#include "TransportCableResult.h"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"



TransportCableResult::TransportCableResult()
{
	rope2.set_Down(true);
}


TransportCableResult::~TransportCableResult()
{
}

void TransportCableResult::prepareSupport()
{
	supports.clear();
	rope1.prepareSupport(supports);
	rope2.prepareSupport(supports);
}

void TransportCableResult::prepareIneqConstraint(int nb_node)
{
	g.resize(nb_node, 1);
	act.resize(nb_node, -1);
	eye_c.resize(nb_node, 0.);
	G.resize(nb_node);	
	for (auto &gg : G) {
		gg.resize(nb_node);
	}
	T.resize(nb_node);
	for (auto &tt : T) {
		tt.resize(nb_node);
	}
	blocked.resize(nb_node);
	blocked_value.resize(nb_node);
	
	KT.resize(nb_node*3);
	for (auto &mm : KT) {
		mm.resize(nb_node*3);
	}
	fi.resize(nb_node*3);
}

int TransportCableResult::exportTC(const std::string & a_fileName, ojson & a_output, const std::string & a_option)
{
	int res = EXIT_SUCCESS;
	try
	{
		res = to_json(a_output, a_option);
	}
	catch (const json::exception& ex) {
		// "Error exporting model " << ex.what());
		throw ex;
	}
	if (a_fileName != "") {
		res = EXIT_FAILURE;
		std::ofstream outputFile(a_fileName);
		if (outputFile.is_open()) {
			outputFile << a_output.dump(2);
			res = EXIT_SUCCESS;
		}
	}
	return res;
}

int TransportCableResult::to_json(ojson & j, const std::string & a_option)
{
	if (a_option == "fem") {		
		j["g"] = g;
		j["eye_c"] = eye_c;
		/*j["G"] = G;		!!trop gros!!
		j["T"] = T;*/

		j["blocked"] = blocked;
		j["blocked_value"] = blocked_value;
	}
	else if (a_option == "cableDS") {
		// M, b, fi, KT
		j["elem_length"] = elem_length;
		//j["b"] = b;
		j["fi"] = fi;
	}
	else if (a_option == "ropeway") {
		rope1.to_json(j["rope1"]);
		rope2.to_json(j["rope2"]);
	}
	else {
		j["supports"] = supports;
		j["pulleys"] = ojson::array();
		j["pulleys"].push_back(puller12);
		j["pulleys"].push_back(puller21);

		j["q"] = q;
		j["punct"] = punct;
		j["act"] = act;

	}
	return EXIT_SUCCESS;
}
