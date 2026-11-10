#pragma once
#include <cuda_runtime.h>
#include "Types.h"
#include "GridInfo.h"

// siconos: no CUDA_MAX_NEIGHBORS cap. kComputeCounts counts every neighbor
// and the lists are allocated from these counts: a cap only dropped the
// neighbors visited last, a subset depending on the order of the points in
// their cells (atomicAdd), hence missed contacts from one run to another.
#define CUDA_META_GRID_GROUP_SIZE 8
#define CUDA_META_GRID_BLOCK_SIZE (CUDA_META_GRID_GROUP_SIZE*CUDA_META_GRID_GROUP_SIZE*CUDA_META_GRID_GROUP_SIZE)

typedef unsigned int uint;
typedef unsigned char byte;
using namespace cuNSearch;

__global__ void kComputeMinMax(
	const Real3 *particles,
	uint particleCount,
	float searchRadius,
	Int3 *minCell,
	Int3 *maxCell
);

__global__ void kInsertParticles(
	const GridInfo GridInfo,
	const Real3 *particles,
	uint *particleCellIndices,
	uint *cellParticleCounts,
	uint *sortIndices
);

__global__ void kInsertParticles_Morton(
	const GridInfo GridInfo,
	const Real3 *particles,
	uint *particleCellIndices,
	uint *cellParticleCounts,
	uint *sortIndices
);

__global__ void kCountingSortIndices(
	const GridInfo GridInfo,
	const uint *particleCellIndices,
	const uint *cellOffsets,
	const uint *sortIndicesSrc,
	uint *sortIndicesDest
);

__global__ void kComputeCounts(
	const Real3 *queryPoints,
	const uint queryPointCount,

	const GridInfo GridInfo,
	const Real3 *particles,
	const uint *cellOffsets,
	const uint *cellParticleCounts,
	uint *neighborCounts,
	const uint *reversedSortIndices,
	const bool selfQuery
);

__global__ void kNeighborhoodQueryWithCounts(
	const Real3 *queryPoints,
	const uint queryPointCount,

	const GridInfo GridInfo,
	const Real3 *particles,
	const uint *cellOffsets,
	const uint *cellParticleCounts,
	const uint *neighborWriteOffsets,
	uint *neighbors,
	const uint *reversedSortIndices,
	const bool selfQuery);
