#include "PointSet.h"
#include "PointSetImplementation.h"
#include "cuda_helper.h"
#include "NotImplementedException.h"

namespace cuNSearch
{
	PointSet::PointSet(PointSet const& other)
	{
		this->m_dynamic = other.m_dynamic;
		this->m_x = other.m_x;
		this->m_n = other.m_n;
		this->m_user_data = other.m_user_data;
		this->sortIndices = other.sortIndices;
		this->neighbors = other.neighbors;

		PointSetImplementation *ptr = other.impl.get();
		impl = make_unique<PointSetImplementation>(PointSetImplementation(*ptr));
	}

	PointSet::~PointSet()
	{

	}

	PointSet::PointSet(Real const* x, std::size_t n, bool dynamic, void *user_data)
		: m_x(x), m_n(n), m_dynamic(dynamic), m_user_data(user_data)
	{
		impl = make_unique<PointSetImplementation>(n, (Real3*)x);
	}

	void PointSet::resize(Real const* x, std::size_t n)
	{
		m_x = x;
		m_n = n;

		impl->resize(n, (Real3*)x);
	}
};
