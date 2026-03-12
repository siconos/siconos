#pragma once
#include <cuda_runtime.h>
#include "Types.h"

__host__ __device__ inline uint CellIndicesToLinearIndex(
	UInt3 &cellDimensions, 
	UInt3 &xyz
)
{
	return xyz.z * cellDimensions.y * cellDimensions.x + xyz.y * cellDimensions.x + xyz.x;
}

__host__ __device__ inline uint CellIndicesToLinearIndex(
	const UInt3&cellDimensions, 
	Int3 &xyz
)
{
	return xyz.z * cellDimensions.y * cellDimensions.x + xyz.y * cellDimensions.x + xyz.x;
}

__host__ __device__ inline void LinearCellIndexTo3DIndices(
	const UInt3 &cellDimensions,
	const uint linearIndex,
	UInt3 &xyz
)
{
	xyz.z = linearIndex / (cellDimensions.y * cellDimensions.x);
	xyz.y = (linearIndex % (cellDimensions.y * cellDimensions.x)) / (cellDimensions.x);
	xyz.x = (linearIndex % (cellDimensions.y * cellDimensions.x)) % cellDimensions.x;
}

__host__ __device__ inline UInt3 LinearCellIndexTo3DIndices(
	const UInt3 &cellDimensions,
	const uint linearIndex)
{
	UInt3 xyz;
	xyz.z = linearIndex / (cellDimensions.y * cellDimensions.x);
	xyz.y = (linearIndex % (cellDimensions.y * cellDimensions.x)) / (cellDimensions.x);
	xyz.x = (linearIndex % (cellDimensions.y * cellDimensions.x)) % cellDimensions.x;
	return xyz;
}

__host__ __device__ inline Int3 LinearCellIndexTo3DIndicesInt3(
	const UInt3 &cellDimensions, 
	const uint &linearIndex)
{
	Int3 xyz;
	xyz.z = linearIndex / (cellDimensions.y * cellDimensions.x);
	xyz.y = (linearIndex % (cellDimensions.y * cellDimensions.x)) / (cellDimensions.x);
	xyz.x = (linearIndex % (cellDimensions.y * cellDimensions.x)) % cellDimensions.x;
	return xyz;
}