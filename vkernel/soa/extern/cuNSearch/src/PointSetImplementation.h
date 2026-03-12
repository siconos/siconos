#pragma once
#include "Types.h"
#include "GridInfo.h"
#include <thrust/device_vector.h>
#include "cuda_helper.h"

namespace cuNSearch
{
	class NeighborhoodSearch;
	class cuNSearchDeviceData;

	class PointSetImplementation
	{
	public:
		PointSetImplementation(size_t particleCount, Real3 *particles);

		PointSetImplementation(PointSetImplementation const& other) = default;
		PointSetImplementation& operator=(PointSetImplementation const& other);
		~PointSetImplementation() { }

		void resize(size_t particleCount, Real3 *particles)
		{
			m_ParticleCount = particleCount;
			m_Particles = particles;

			uint threadStarts = 0;
			CudaHelper::GetThreadBlocks(static_cast<unsigned int>(particleCount), ThreadsPerBlock, BlockStartsForParticles, threadStarts);

			copyToDevice();
		}

		void copyToDevice();

	private:
		friend NeighborhoodSearch;
		friend cuNSearchDeviceData;

		// Min Max of all particles
		Real3 Min, Max;

		size_t m_ParticleCount;

		// Pointer to the host particle data
		Real3 *m_Particles;

		// Number of thread blocks that must be started to start a thread per particle
		int ThreadsPerBlock;
		uint BlockStartsForParticles;

		// All device data for the a point set to perform query operations.
		GridInfo gridInfo;
		thrust::device_vector<Real3> d_Particles;
		thrust::device_vector<uint> d_ParticleCellIndices;

		thrust::device_vector<uint> d_CellOffsets;
		thrust::device_vector<uint> d_CellParticleCounts;
		thrust::device_vector<uint> d_SortIndices;
		thrust::device_vector<uint> d_ReversedSortIndices;

		void prepareInternalDataStructures(GridInfo &gridInfo, size_t numberOfCells);
	};
};