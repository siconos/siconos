#pragma once
#include "Types.h"
#include "PointSet.h"
#include <thrust/device_vector.h>

namespace cuNSearch
{
	class cuNSearchDeviceData
	{
	public:
		cuNSearchDeviceData(Real searchRadius)
		{
			m_SearchRadius = searchRadius;
		}

		void setSearchRadius(Real searchRadius)
		{
			m_SearchRadius = searchRadius;
		}

		/** Compute min max for the given point set
		*/
		void computeMinMax(PointSet &pointSet);


		/** Constructs the uniform grid for the given point set a updates all cell information to allow queries on this point set
		*/
		void computeCellInformation(PointSet &pointSet);


		/** Queries the neighbors in the given point set for all particles in the query point set.
		*/
		void computeNeighborhood(PointSet &queryPointSet, PointSet &pointSet, uint neighborListEntry);

	private:
		Real m_SearchRadius;

		thrust::device_vector<Int3> d_MinMax;
		thrust::device_vector<uint> d_TempSortIndices;

		//Device neighbor buffers (only temporary used: after the computation the data is copied to the host)
		thrust::device_vector<uint> d_Neighbors;
		thrust::device_vector<uint> d_NeighborCounts;
		thrust::device_vector<uint> d_NeighborWriteOffsets;
	};
};