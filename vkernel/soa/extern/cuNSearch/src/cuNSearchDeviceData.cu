#include "cuNSearchDeviceData.h"

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

#include "PointSetImplementation.h"
#include "GridInfo.h"
#include "cuda_helper.h"
#include "cuNSearchKernels.cuh"

namespace cuNSearch
{
	void cuNSearchDeviceData::computeMinMax(PointSet &pointSet)
	{
		if (pointSet.n_points() == 0)
			return;
			
		auto pointSetImpl = pointSet.impl.get();

		Int3 data[2];
		data[0] = Int3(std::numeric_limits<int>().max(), std::numeric_limits<int>().max(), std::numeric_limits<int>().max());
		data[1] = Int3(std::numeric_limits<int>().min(), std::numeric_limits<int>().min(), std::numeric_limits<int>().min());
		d_MinMax.resize(2);
		CudaHelper::MemcpyHostToDevice(data, CudaHelper::GetPointer(d_MinMax), 2);

		kComputeMinMax << <pointSetImpl->BlockStartsForParticles, pointSetImpl->ThreadsPerBlock >> > (
			(Real3*)CudaHelper::GetPointer(pointSetImpl->d_Particles),
			static_cast<unsigned int>(pointSet.n_points()),
			m_SearchRadius,
			CudaHelper::GetPointer(d_MinMax), 
			CudaHelper::GetPointer(d_MinMax) + 1
			);
		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_MinMax), data, 2);
		Int3 minCell = data[0];
		Int3 maxCell = data[1];

		pointSetImpl->Min.x = minCell.x * m_SearchRadius;
		pointSetImpl->Min.y = minCell.y * m_SearchRadius;
		pointSetImpl->Min.z = minCell.z * m_SearchRadius;

		pointSetImpl->Max.x = maxCell.x * m_SearchRadius;
		pointSetImpl->Max.y = maxCell.y * m_SearchRadius;
		pointSetImpl->Max.z = maxCell.z * m_SearchRadius;

		//CPU implementation of min max computation
		//Real3 cpuMin, cpuMax;
		//cpuMin = make_Real3(std::numeric_limits<Real>().max());
		//cpuMax = make_Real3(std::numeric_limits<Real>().min());

		//Real3 *points = (Real3 *)pointSet.m_x;
		//for (size_t i = 0; i < pointSet.n_points(); i++)
		//{
		//	cpuMin.x = std::min(cpuMin.x, points[i].x);
		//	cpuMin.y = std::min(cpuMin.y, points[i].y);
		//	cpuMin.z = std::min(cpuMin.z, points[i].z);

		//	cpuMax.x = std::max(cpuMax.x, points[i].x);
		//	cpuMax.y = std::max(cpuMax.y, points[i].y);
		//	cpuMax.z = std::max(cpuMax.z, points[i].z);
		//}
	}

	void cuNSearchDeviceData::computeCellInformation(PointSet &pointSet)
	{
		if (pointSet.n_points() == 0)
			return;
	
		auto pointSetImpl = pointSet.impl.get();
		Real3 sceneMin = pointSetImpl->Min;
		Real3 sceneMax = pointSetImpl->Max;

		GridInfo gridInfo;
		gridInfo.ParticleCount = static_cast<uint>(pointSet.n_points());
		gridInfo.SquaredSearchRadius = m_SearchRadius * m_SearchRadius;
		gridInfo.GridMin = sceneMin;

		Real cellSize = m_SearchRadius;
		Real3 gridSize = sceneMax - sceneMin;
		gridInfo.GridDimension.x = static_cast<unsigned int>(ceil(gridSize.x / cellSize));
		gridInfo.GridDimension.y = static_cast<unsigned int>(ceil(gridSize.y / cellSize));
		gridInfo.GridDimension.z = static_cast<unsigned int>(ceil(gridSize.z / cellSize));

		//Increase grid by 2 cells in each direciton (+4 in each dimension) to skip bounds checks in the kernel
		gridInfo.GridDimension.x += 4;
		gridInfo.GridDimension.y += 4;
		gridInfo.GridDimension.z += 4;
		gridInfo.GridMin -= Real3(cellSize, cellSize, cellSize) * (Real)2;

		//One meta grid cell contains 8x8x8 grild cells. (512)
		gridInfo.MetaGridDimension.x = static_cast<unsigned int>(ceil(gridInfo.GridDimension.x / (float)CUDA_META_GRID_GROUP_SIZE));
		gridInfo.MetaGridDimension.y = static_cast<unsigned int>(ceil(gridInfo.GridDimension.y / (float)CUDA_META_GRID_GROUP_SIZE));
		gridInfo.MetaGridDimension.z = static_cast<unsigned int>(ceil(gridInfo.GridDimension.z / (float)CUDA_META_GRID_GROUP_SIZE));

		// Adjust grid size to multiple of cell size
		gridSize.x = gridInfo.GridDimension.x * cellSize;
		gridSize.y = gridInfo.GridDimension.y * cellSize;
		gridSize.z = gridInfo.GridDimension.z * cellSize;

		gridInfo.GridDelta.x = gridInfo.GridDimension.x / gridSize.x;
		gridInfo.GridDelta.y = gridInfo.GridDimension.y / gridSize.y;
		gridInfo.GridDelta.z = gridInfo.GridDimension.z / gridSize.z;

		d_TempSortIndices.resize(gridInfo.ParticleCount);

		uint numberOfCells = (gridInfo.MetaGridDimension.x * gridInfo.MetaGridDimension.y * gridInfo.MetaGridDimension.z) * CUDA_META_GRID_BLOCK_SIZE;
		pointSet.impl->prepareInternalDataStructures(gridInfo, numberOfCells);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		cudaMemset(CudaHelper::GetPointer(pointSetImpl->d_CellParticleCounts), 0, CudaHelper::GetSizeInBytes(pointSetImpl->d_CellParticleCounts));

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		kInsertParticles_Morton << <pointSetImpl->BlockStartsForParticles, pointSetImpl->ThreadsPerBlock >> > (
			gridInfo,
			(Real3*)CudaHelper::GetPointer(pointSetImpl->d_Particles),
			CudaHelper::GetPointer(pointSetImpl->d_ParticleCellIndices),
			CudaHelper::GetPointer(pointSetImpl->d_CellParticleCounts),
			CudaHelper::GetPointer(d_TempSortIndices)
			);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		thrust::exclusive_scan(
			pointSetImpl->d_CellParticleCounts.begin(),
			pointSetImpl->d_CellParticleCounts.end(),
			pointSetImpl->d_CellOffsets.begin());
		CudaHelper::DeviceSynchronize();

		kCountingSortIndices << <pointSetImpl->BlockStartsForParticles, pointSetImpl->ThreadsPerBlock >> > (
			gridInfo,
			CudaHelper::GetPointer(pointSetImpl->d_ParticleCellIndices),
			CudaHelper::GetPointer(pointSetImpl->d_CellOffsets),
			CudaHelper::GetPointer(d_TempSortIndices),
			CudaHelper::GetPointer(pointSetImpl->d_SortIndices)
			);

		CudaHelper::DeviceSynchronize();

		auto &tempSequence = d_TempSortIndices;
		thrust::sequence(tempSequence.begin(), tempSequence.end());

		thrust::gather(
			pointSetImpl->d_SortIndices.begin(),
			pointSetImpl->d_SortIndices.end(),
			tempSequence.begin(),
			pointSetImpl->d_ReversedSortIndices.begin());

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		pointSet.sortIndices.resize(pointSetImpl->d_SortIndices.size());
		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(pointSetImpl->d_SortIndices), pointSet.sortIndices.data(), pointSetImpl->d_SortIndices.size());
	}

	void cuNSearchDeviceData::computeNeighborhood(PointSet &queryPointSet, PointSet &pointSet, uint neighborListEntry)
	{
		if (queryPointSet.n_points() == 0)
			return;
	
		auto queryPointSetImpl = queryPointSet.impl.get();
		auto pointSetImpl = pointSet.impl.get();

		uint particleCount = static_cast<uint>(queryPointSet.n_points());

		USE_TIMING(Timing::startTiming("Execute kNeighborCount"));
		d_NeighborCounts.resize(particleCount);

		kComputeCounts << <queryPointSetImpl->BlockStartsForParticles, queryPointSetImpl->ThreadsPerBlock >> > (
			(Real3*)CudaHelper::GetPointer(queryPointSetImpl->d_Particles),
			static_cast<unsigned int>(queryPointSet.n_points()),

			pointSetImpl->gridInfo,
			(Real3*)CudaHelper::GetPointer(pointSetImpl->d_Particles),
			CudaHelper::GetPointer(pointSetImpl->d_CellOffsets),
			CudaHelper::GetPointer(pointSetImpl->d_CellParticleCounts),

			CudaHelper::GetPointer(d_NeighborCounts),
			CudaHelper::GetPointer(pointSetImpl->d_ReversedSortIndices)
			);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();

		USE_TIMING(Timing::stopTiming(PRINT_STATS));
		USE_TIMING(Timing::startTiming("Execute exclusive_scan over counts"));

		d_NeighborWriteOffsets.resize(particleCount);

		//Prefix sum over neighbor counts
		thrust::exclusive_scan(
			d_NeighborCounts.begin(),
			d_NeighborCounts.end(),
			d_NeighborWriteOffsets.begin());

		CudaHelper::DeviceSynchronize();

		//Compute total amount of neighbors
		uint lastOffset = 0;
		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_NeighborWriteOffsets) + particleCount - 1, &lastOffset, 1);
		uint lastParticleNeighborCount = 0;
		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_NeighborCounts) + particleCount - 1, &lastParticleNeighborCount, 1);
		uint totalNeighborCount = lastOffset + lastParticleNeighborCount;
		d_Neighbors.resize(totalNeighborCount);

		CudaHelper::DeviceSynchronize();

		USE_TIMING(Timing::stopTiming(PRINT_STATS));
		USE_TIMING(Timing::startTiming("Execute kNeighborhoodQueryWithCounts"));

		kNeighborhoodQueryWithCounts << <queryPointSetImpl->BlockStartsForParticles, queryPointSetImpl->ThreadsPerBlock >> > (
			(Real3*)CudaHelper::GetPointer(queryPointSetImpl->d_Particles),
			static_cast<unsigned int>(queryPointSet.n_points()),

			pointSetImpl->gridInfo,
			(Real3*)CudaHelper::GetPointer(pointSetImpl->d_Particles),
			CudaHelper::GetPointer(pointSetImpl->d_CellOffsets),
			CudaHelper::GetPointer(pointSetImpl->d_CellParticleCounts),

			CudaHelper::GetPointer(d_NeighborWriteOffsets),
			CudaHelper::GetPointer(d_Neighbors),
			CudaHelper::GetPointer(pointSetImpl->d_ReversedSortIndices)
			);

		CudaHelper::CheckLastError();
		CudaHelper::DeviceSynchronize();
		USE_TIMING(Timing::stopTiming(PRINT_STATS));

		//Copy data to host
		USE_TIMING(Timing::startTiming("Neighbor copy from device to host - resize"));

		auto &neighborSet = queryPointSet.neighbors[neighborListEntry];

		if (neighborSet.NeighborCountAllocationSize < totalNeighborCount)
		{
			if (neighborSet.NeighborCountAllocationSize != 0)
			{
				cudaFreeHost(neighborSet.Neighbors);
			}

			neighborSet.NeighborCountAllocationSize = static_cast<unsigned int>(totalNeighborCount * 1.5);
			cudaMallocHost(&neighborSet.Neighbors, sizeof(uint) * neighborSet.NeighborCountAllocationSize);
		}
		if (neighborSet.ParticleCountAllocationSize < particleCount)
		{
			if (neighborSet.ParticleCountAllocationSize != 0)
			{
				cudaFreeHost(neighborSet.Offsets);
				cudaFreeHost(neighborSet.Counts);
			}

			neighborSet.ParticleCountAllocationSize = static_cast<unsigned int>(particleCount * 1.5);
			cudaMallocHost(&neighborSet.Offsets, sizeof(uint) * neighborSet.ParticleCountAllocationSize);
			cudaMallocHost(&neighborSet.Counts, sizeof(uint) * neighborSet.ParticleCountAllocationSize);
		}

		USE_TIMING(Timing::stopTiming(PRINT_STATS));
		USE_TIMING(Timing::startTiming("Neighbor copy from device to host - MemcpyDeviceToHost"));

		if (PRINT_STATS)
		{
			int bytesToCopy = totalNeighborCount * 4 + particleCount * 2 * 4;
			printf("Total neighbors: %d \n", totalNeighborCount);
			printf("Average neighbors: %d \n", totalNeighborCount / particleCount);
			printf("Expected amount: %f MB \n", bytesToCopy / (1024.0f * 1024.0f));
		}

		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_Neighbors), neighborSet.Neighbors, totalNeighborCount);
		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_NeighborCounts), neighborSet.Counts, particleCount);
		CudaHelper::MemcpyDeviceToHost(CudaHelper::GetPointer(d_NeighborWriteOffsets), neighborSet.Offsets, particleCount);

		USE_TIMING(Timing::stopTiming(PRINT_STATS));
	}
}