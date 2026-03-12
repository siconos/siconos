
#include "cuNSearch.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <chrono>
#include <thread>

#include <cuda_runtime.h>

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/gather.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/fill.h>

#ifdef DEBUG
#define PRINT_STATS true
#define USE_TIMING(x) x;
#else
#define PRINT_STATS false
#define USE_TIMING(x)
#endif
#include "Timing.h"

#include "NotImplementedException.h"
#include "cuNSearchDeviceData.h"
#include "PointSetImplementation.h"
#include "cuNSearchKernels.cuh"
#include "cuda_helper.h"


namespace cuNSearch
{
	NeighborhoodSearch::NeighborhoodSearch(Real searchRadius)
	{
		deviceData = make_unique<cuNSearchDeviceData>(searchRadius);
		set_radius(searchRadius);
	}

	NeighborhoodSearch::~NeighborhoodSearch()
	{

	}

	unsigned int NeighborhoodSearch::add_point_set(Real const* x, std::size_t n, bool is_dynamic,
		bool search_neighbors, bool find_neighbors, void *user_data)
	{
		auto index = pointSets.size();
		pointSets.push_back(PointSet(x, n, is_dynamic, user_data));
		m_activation_table.add_point_set(search_neighbors, find_neighbors);

		for (auto &pointSet : pointSets)
		{
			pointSet.neighbors.resize(pointSets.size());
		}

		return static_cast<unsigned int>(index);
	}


	void NeighborhoodSearch::set_radius(Real r)
	{
		this->searchRadius = r;
		deviceData->setSearchRadius(r);
		isInitialized = false;
	}

	void NeighborhoodSearch::z_sort()
	{
		//Do nothing as the sort step is part of the main implementation
	}

	void
		NeighborhoodSearch::resize_point_set(unsigned int index, Real const* x, std::size_t size)
	{
		pointSets[index].resize(x, size);
	}

	void
		NeighborhoodSearch::update_activation_table()
	{
		//Update neighborhood search data structures after changing the activation table.
		//If general find_neighbors() function is called there is no requirement to manually update the point sets.
	}

	void
		NeighborhoodSearch::updatePointSet(PointSet &pointSet)
	{
		USE_TIMING(Timing::startTiming("Update point sets - copyParticleData"));
		pointSet.impl->copyToDevice();
		USE_TIMING(Timing::stopTiming(PRINT_STATS));

		USE_TIMING(Timing::startTiming("Update point sets - computeMinMax"));
		deviceData->computeMinMax(pointSet);
		USE_TIMING(Timing::stopTiming(PRINT_STATS));

		USE_TIMING(Timing::startTiming("Update point sets - computeCellInformation"));
		deviceData->computeCellInformation(pointSet);

		USE_TIMING(Timing::stopTiming(PRINT_STATS));
	}

	void
		NeighborhoodSearch::find_neighbors(bool points_changed_)
	{
		if (points_changed_ || !isInitialized)
		{
			for (auto &pointSet : pointSets)
			{
				if (!isInitialized || pointSet.is_dynamic())
				{
					updatePointSet(pointSet);
				}
			}
		}
		isInitialized = true;

		for (unsigned int i = 0; i < pointSets.size(); i++)
		{
			for (unsigned int j = 0; j < pointSets.size(); j++)
			{
				if (m_activation_table.is_active(i, j))
				{
					auto &queryPointSet = pointSets[i];
					auto &pointSet = pointSets[j];
					deviceData->computeNeighborhood(queryPointSet, pointSet, j);
				}
			}
		}
	}

	void
		NeighborhoodSearch::find_neighbors(unsigned int point_set_id, unsigned int point_index, std::vector<std::vector<unsigned int>> &neighbors)
	{
		throw new NotImplementedException("NeighborhoodSearch::find_neighbors()");
	}

	void
		NeighborhoodSearch::update_point_sets()
	{
		for (unsigned int i = 0; i < pointSets.size(); i++)
		{
			update_point_set(i);
		}
	}

	void
		NeighborhoodSearch::update_point_set(int i)
	{
		updatePointSet(pointSets[i]);
	}
}

