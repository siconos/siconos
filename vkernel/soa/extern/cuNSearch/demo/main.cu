#include "cuNSearch.h"
#include "Timing.h"

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <random>

using namespace cuNSearch;

using Real3 = std::array<Real, 3>;
std::vector<Real3> positions;

inline Real3 operator-(const Real3 & left, const Real3 & right)
{
	return Real3{ left[0] - right[0], left[1] - right[1], left[2] - right[2] };
}

std::size_t const N = 120;
Real const r_omega = static_cast<Real>(0.15);
Real const r_omega2 = r_omega * r_omega;
Real const radius = static_cast<Real>(2.0) * (static_cast<Real>(2.0) * r_omega / static_cast<Real>(N - 1));

void testCuNSearch()
{
	//Generate test data
	Real min_x = std::numeric_limits<Real>::max();
	Real max_x = std::numeric_limits<Real>::min();
	positions.reserve(N * N * N);
	for (unsigned int i = 0; i < N; ++i)
	{
		for (unsigned int j = 0; j < N; ++j)
		{
			for (unsigned int k = 0; k < N; ++k)
			{
				std::array<Real, 3> x = { {
						r_omega * static_cast<Real>(2.0 * static_cast<double>(i) / static_cast<double>(N - 1) - 1.0),
						r_omega * static_cast<Real>(2.0 * static_cast<double>(j) / static_cast<double>(N - 1) - 1.0),
						r_omega * static_cast<Real>(2.0 * static_cast<double>(k) / static_cast<double>(N - 1) - 1.0) } };

				Real l2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
				if (l2 < r_omega2)
				{
					x[0] += static_cast<Real>(0.35);
					x[1] += static_cast<Real>(0.35);
					x[2] += static_cast<Real>(0.35);
					positions.push_back(x);
					if (min_x > x[0])
					{
						min_x = x[0];
					}
					if (max_x < x[0])
					{
						max_x = x[0];
					}
				}
			}
		}
	}
	std::random_shuffle(positions.begin(), positions.end());
	printf("Number of particles: %d \n", static_cast<int>(positions.size()));

	//Create neighborhood search instance
	NeighborhoodSearch nsearch(radius);

	//Add point set from the test data
	auto pointSetIndex = nsearch.add_point_set(positions.front().data(), positions.size(), true, true);

	for (size_t i = 0; i < 5; i++)
	{
		if (i != 0)
		{
			nsearch.z_sort();
			nsearch.point_set(pointSetIndex).sort_field((Real3*)nsearch.point_set(pointSetIndex).GetPoints());
		}

		Timing::reset();
		nsearch.find_neighbors();
		Timing::printAverageTimes();
	}

	//Neighborhood search result test
	auto &pointSet = nsearch.point_set(0);
	auto points = pointSet.GetPoints();

	std::cout << "Validate results" << std::endl;
	for (unsigned int i = 0; i < pointSet.n_points(); i++)
	{
		Real3 point = ((Real3*)points)[i];
		auto count = pointSet.n_neighbors(0, i);
		for (unsigned int j = 0; j < count; j++)
		{
			auto neighbor = pointSet.neighbor(0, i, j);
			auto diff = point - ((Real3*)points)[neighbor];
			float squaredLength = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
			float distance = sqrt(squaredLength);

			if (distance > radius)
			{
				throw std::runtime_error("Not a neighbor");
			}
		}
	}
}

int main(int argc, char* argv[])
{
#ifdef DEBUG
	std::cout << "Debug Build:" << std::endl;

	if(sizeof(Real) == 4)
		std::cout << "Real = float" << std::endl;
	else if (sizeof(Real) == 8)
		std::cout << "Real = double" << std::endl;
#endif

	testCuNSearch();
	std::cout << "Finished Testing" << std::endl;
	getchar();
}
