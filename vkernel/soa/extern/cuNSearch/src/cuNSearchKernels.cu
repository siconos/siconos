#include "cuNSearchKernels.cuh"
#include "Types.h"
#include "helper_mortonCode.h"
#include "helper_linearIndex.h"

#define INT16_RANGE 32767
#define UPDATE_REF_OFFSET -32768

#pragma region HelperMethods
inline __device__ uint ToCellIndex_MortonMetaGrid(const GridInfo &GridInfo, Int3 gridCell)
{
	Int3 metaGridCell = Int3(
		gridCell.x / CUDA_META_GRID_GROUP_SIZE,
		gridCell.y / CUDA_META_GRID_GROUP_SIZE,
		gridCell.z / CUDA_META_GRID_GROUP_SIZE);
	gridCell.x %= CUDA_META_GRID_GROUP_SIZE;
	gridCell.y %= CUDA_META_GRID_GROUP_SIZE;
	gridCell.z %= CUDA_META_GRID_GROUP_SIZE;
	uint metaGridIndex = CellIndicesToLinearIndex(GridInfo.MetaGridDimension, metaGridCell);
	return metaGridIndex * CUDA_META_GRID_BLOCK_SIZE + MortonCode3(gridCell.x, gridCell.y, gridCell.z);
}

inline __device__ Int3 ToGridCell_MortonMetaGrid(const GridInfo &GridInfo, uint &cellIndex)
{
	uint metaGridIndex = cellIndex / CUDA_META_GRID_BLOCK_SIZE;
	Int3 gridCell = MortonCodeToIndexInt3(cellIndex % CUDA_META_GRID_BLOCK_SIZE) + LinearCellIndexTo3DIndicesInt3(GridInfo.MetaGridDimension, metaGridIndex) * CUDA_META_GRID_GROUP_SIZE;
	return gridCell;
}

inline __device__ uint ToCellIndex(const GridInfo &GridInfo, Int3 gridCell)
{
	return CellIndicesToLinearIndex(GridInfo.GridDimension, gridCell);
}

inline __device__ Int3 ToGridCell(const GridInfo &GridInfo, const uint &cellIndex)
{
	return LinearCellIndexTo3DIndicesInt3(GridInfo.GridDimension, cellIndex);
}
#pragma endregion HelperMethods

__global__ void kComputeMinMax(
	const Real3 *particles,
	uint particleCount,
	float searchRadius,
	Int3 *minCell,
	Int3 *maxCell
)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= particleCount) return;
	const Real3 particle = particles[particleIndex];

	Int3 cell;
	cell.x = (int)floor(particle.x / searchRadius);
	cell.y = (int)floor(particle.y / searchRadius);
	cell.z = (int)floor(particle.z / searchRadius);

	atomicMin(&(minCell->x), cell.x);
	atomicMin(&(minCell->y), cell.y);
	atomicMin(&(minCell->z), cell.z);

	atomicMax(&(maxCell->x), cell.x);
	atomicMax(&(maxCell->y), cell.y);
	atomicMax(&(maxCell->z), cell.z);

	//printf("%d %d %d Min: %d %d %d Max: %d %d %d \n", cell.x, cell.y, cell.z, minCell->x, minCell->y, minCell->z, maxCell->x, maxCell->y, maxCell->z);
}


#pragma region kInsertParticles
__global__ void kInsertParticles(
	const GridInfo GridInfo,
	const Real3 *particles,
	uint *particleCellIndices,
	uint *cellParticleCounts,
	uint *sortIndices
)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= GridInfo.ParticleCount) return;

	Real3 gridCellF = (particles[particleIndex] - GridInfo.GridMin) * GridInfo.GridDelta;
	Int3 gridCell = Int3(int(gridCellF.x), int(gridCellF.y), int(gridCellF.z));
	uint cellIndex = ToCellIndex(GridInfo, gridCell);
	particleCellIndices[particleIndex] = cellIndex;
	sortIndices[particleIndex] = atomicAdd(&cellParticleCounts[cellIndex], 1);
}

__global__ void kInsertParticles_Morton(
	const GridInfo GridInfo,
	const Real3 *particles,
	uint *particleCellIndices,
	uint *cellParticleCounts,
	uint *sortIndices
)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= GridInfo.ParticleCount) return;

	Real3 gridCellF = (particles[particleIndex] - GridInfo.GridMin) * GridInfo.GridDelta;
	Int3 gridCell = Int3(int(gridCellF.x), int(gridCellF.y), int(gridCellF.z));
	uint cellIndex = ToCellIndex_MortonMetaGrid(GridInfo, gridCell);
	particleCellIndices[particleIndex] = cellIndex;
	sortIndices[particleIndex] = atomicAdd(&cellParticleCounts[cellIndex], 1);
}
#pragma endregion kInsertParticles

__global__ void kCountingSortIndices(
	const GridInfo GridInfo,
	const uint *particleCellIndices,
	const uint *cellOffsets,
	const uint *sortIndicesSrc,
	uint *sortIndicesDest
)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= GridInfo.ParticleCount) return;

	uint gridCellIndex = particleCellIndices[particleIndex];

	uint sortIndex = sortIndicesSrc[particleIndex] + cellOffsets[gridCellIndex];
	sortIndicesDest[sortIndex] = particleIndex;
}

__global__ void kComputeCounts(
	const Real3 *queryPoints,
	const uint queryPointCount,

	const GridInfo GridInfo,
	const Real3 *particles,
	const uint *cellOffsets,
	const uint *cellParticleCounts,
	uint *neighborCounts,
	const uint *reversedSortIndices
)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= queryPointCount) return;
	const Real3 particle = queryPoints[particleIndex];
	Real3 gridCellF = (particle - GridInfo.GridMin) * GridInfo.GridDelta;

	Int3 coord = Int3(int(floor(gridCellF.x)), int(floor(gridCellF.y)), int(floor(gridCellF.z)));

	uint neighborCount = 0;
	for (int z = -1; z < 2; z++)
		for (int y = -1; y < 2; y++)
			for (int x = -1; x < 2; x++)
			{
				Int3 finalCoord = coord + Int3(x, y, z);

				if (finalCoord.x < 0 || finalCoord.y < 0 || finalCoord.z < 0
					|| finalCoord.x >= GridInfo.GridDimension.x || finalCoord.y >= GridInfo.GridDimension.y || finalCoord.z >= GridInfo.GridDimension.z)
					continue;

				uint neighborCellIndex = ToCellIndex_MortonMetaGrid(GridInfo, finalCoord);
				uint neighborCellCount = cellParticleCounts[neighborCellIndex];
				uint neighborCellStart = cellOffsets[neighborCellIndex];

				for (uint i = neighborCellStart; i < neighborCellStart + neighborCellCount; i++)
				{
					uint &neighborIndex = i;
					Real3 diff = particles[reversedSortIndices[neighborIndex]] - particle;
					float squaredDistance = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
					if (squaredDistance < GridInfo.SquaredSearchRadius && squaredDistance > 0.0)
					{
						neighborCount++;
					}

					if (neighborCount == CUDA_MAX_NEIGHBORS)
					{
						neighborCounts[particleIndex] = neighborCount;
						return;
					}
				}
			}

	neighborCounts[particleIndex] = neighborCount;
}


__global__ void kNeighborhoodQueryWithCounts(
	const Real3 *queryPoints,
	const uint queryPointCount,

	const GridInfo GridInfo,
	const Real3 *particles,
	const uint *cellOffsets,
	const uint *cellParticleCounts,
	const uint *neighborWriteOffsets,
	uint *neighbors,
	const uint *reversedSortIndices)
{
	uint particleIndex = blockIdx.x * blockDim.x + threadIdx.x;
	if (particleIndex >= queryPointCount) return;
	const Real3 particle = queryPoints[particleIndex];
	Real3 gridCellF = (particle - GridInfo.GridMin) * GridInfo.GridDelta;

	Int3 coord = Int3(int(floor(gridCellF.x)), int(floor(gridCellF.y)), int(floor(gridCellF.z)));

	uint neighborCount = 0;
	const uint writeOffset = neighborWriteOffsets[particleIndex];

	for (int z = -1; z < 2; z++)
		for (int y = -1; y < 2; y++)
			for (int x = -1; x < 2; x++)
			{
				Int3 finalCoord = coord + Int3(x, y, z);

				if (finalCoord.x < 0 || finalCoord.y < 0 || finalCoord.z < 0
					|| finalCoord.x >= GridInfo.GridDimension.x || finalCoord.y >= GridInfo.GridDimension.y || finalCoord.z >= GridInfo.GridDimension.z)
					continue;

				uint neighborCellIndex = ToCellIndex_MortonMetaGrid(GridInfo, finalCoord);
				uint neighborCellCount = cellParticleCounts[neighborCellIndex];
				uint neighborCellStart = cellOffsets[neighborCellIndex];

				for (uint i = neighborCellStart; i < neighborCellStart + neighborCellCount; i++)
				{
					uint &neighborIndex = i;
					Real3 diff = particles[reversedSortIndices[neighborIndex]] - particle;
					float squaredDistance = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
					if (squaredDistance < GridInfo.SquaredSearchRadius && squaredDistance > 0.0)
					{
						neighbors[writeOffset + neighborCount] = reversedSortIndices[neighborIndex];
						neighborCount++;
					}

					if (neighborCount == CUDA_MAX_NEIGHBORS)
					{
						return;
					}
				}
			}
}