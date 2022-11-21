#pragma once
#include "Rope.h"
#include "SiconosVector.hpp"

class Support
{
public:
	Support(const Pile &a_pile);
	virtual ~Support();
	virtual const double &get_radius() const;

	virtual void prepare(const Rope &a_rope);
    virtual void prepare(const Pile &a_start, const Pile &a_end, double T);


	virtual void compute(
		const Point &a_p, 
		double a_tol,
		double &g, 
		Point &G,
		Point &T);

	virtual bool isContact(const Point &a_p, const double &a_tol);
	virtual void getContact(const Point &a_p, SiconosVector &pc2, SiconosVector &normal, SiconosVector &tangent);

	friend void to_json(ojson &j, const Support &p);

protected:
	const Pile &r_pile;
	Point m_p;
};

