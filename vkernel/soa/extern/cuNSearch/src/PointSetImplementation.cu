#include "PointSetImplementation.h"
#include "NotImplementedException.h"

namespace cuNSearch
{
	PointSetImplementation::PointSetImplementation(size_t particleCount, Real3 *particles)
	{
		m_ParticleCount = particleCount;
		m_Particles = particles;

		uint threadStarts = 0;
		ThreadsPerBlock = 64;
		CudaHelper::GetThreadBlocks(static_cast<unsigned int>(particleCount), ThreadsPerBlock, BlockStartsForParticles, threadStarts);

		copyToDevice();
	}

	PointSetImplementation& PointSetImplementation::operator=(PointSetImplementation const& other)
	{
		if (this != &other)
		{
			PointSetImplementation tmp(other);
			std::swap(tmp, *this);
		}
		return *this;
	}

	void PointSetImplementation::prepareInternalDataStructures(GridInfo &gridInfo, size_t numberOfCells)
	{
		this->gridInfo = gridInfo;

		d_ParticleCellIndices.resize(m_ParticleCount);
		d_SortIndices.resize(m_ParticleCount);
		d_ReversedSortIndices.resize(m_ParticleCount);

		d_CellOffsets.resize(numberOfCells);
		d_CellParticleCounts.resize(numberOfCells);
	}

	void PointSetImplementation::copyToDevice()
	{
		d_Particles.resize(m_ParticleCount);
		CudaHelper::MemcpyHostToDevice(m_Particles, CudaHelper::GetPointer(d_Particles), m_ParticleCount);
	}
}