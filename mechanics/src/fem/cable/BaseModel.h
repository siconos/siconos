#pragma once
#include "TransportCable.h"

class BaseModel
{
public:
	BaseModel();
	virtual ~BaseModel();

	void from_json(const json &j, const string &a_header);

protected:
	virtual void from_json(const json &j) = 0;

};

