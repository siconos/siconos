#pragma once
#include "TransportCable.h"
#include "Cable.h"
#include "Carriers.h"
#include "Pile.h"

class TransportCableModel
{
public:
	TransportCableModel();
	virtual ~TransportCableModel();

	int from_file(const std::string &a_fileName);
	int from_json(const json &j);
	bool isLoaded();

	const Cable &get_cable() const;
	const Carriers &get_carriers() const;
	const vector<Pile> &get_piles1() const ;
	const vector<Pile> &get_piles2() const;

private:
	Cable m_cable;
	Carriers m_carriers;
	Pile m_stationUp;	// drive station
	Pile m_stationDown;
	vector<Pile> m_piles; // Roller batteries

	void clear();
	int validate();
	vector<Pile> m_pilesUp;
	vector<Pile> m_pilesDown;
};

