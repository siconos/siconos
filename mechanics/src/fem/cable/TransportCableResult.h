#pragma once
#include "Ropeway.h"
#include "Pulley.h"
#include "Support.h"

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

	std::vector<Point> q;	// positions
	std::vector<Point> R;	// internal forces [x,y,z]-> [H,V,B]
	std::vector<double> TS; // tension
	std::vector<int> contacts; // points en contact (=1)

	int nb_nodes;
	double length;
	double elem_length;
	// à convertir en siconos (vecteur ou matrice)
	vector<double> punct;

	vector<double> g;	
	vector<vector<Point>> G;
	vector<vector<Point>> T;


#ifndef NSICONOS
	std::shared_ptr<class SiconosVector> q0{nullptr};
    std::shared_ptr<class SiconosVector> v0{nullptr};

    std::shared_ptr<class SiconosMatrix> mass{nullptr};
    std::shared_ptr<class SiconosVector> b{nullptr};
#endif

};

