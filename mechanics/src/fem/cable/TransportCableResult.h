#pragma once
#include "Ropeway.h"
#include "Pulley.h"
#include "Support.h"
//#include "SiconosVector.hpp"
//#include "SimpleMatrix.hpp"

class TransportCableResult
{
public:
	TransportCableResult();
	virtual ~TransportCableResult();

	void prepareSupport();
	void prepareIneqConstraint(int nb_nodes);

	int exportTC(const std::string &a_fileName,
		ojson &a_output,
		const std::string &a_option = "all");

	int to_json(ojson &j, const std::string &a_option = "all");

	int puller12idx;
	int puller21idx;
	Ropeway rope1;
	Ropeway rope2;

	vector<std::shared_ptr<Support>> supports;

	vector<Point> q;

	double elem_length;
	// à convertir en siconos (vecteur ou matrice)
	vector<double> punct;

	vector<double> g;
	vector<int> act;
	vector<double> eye_c;
	vector<vector<Point>> G;
	vector<vector<Point>> T;
	vector<Point> blocked;
	vector<Point> blocked_value;
	
	std::shared_ptr<class SiconosVector> q0{nullptr};
    std::shared_ptr<class SiconosVector> v0{nullptr};

    std::shared_ptr<class SiconosMatrix> mass{nullptr};
    std::shared_ptr<class SiconosVector> b{nullptr};



	vector<vector<double>> KT;
	vector<double> fi;
};

