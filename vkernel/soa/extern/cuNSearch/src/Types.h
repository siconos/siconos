#pragma once
#include "Common.h"

#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "="  VALUE(var)
#pragma message(VAR_NAME_VALUE(CUDA_VERSION))

#include <cuda.h>
#include <cuda_runtime.h>

#define CUDA_METHOD __host__ __device__

namespace cuNSearch
{
	#ifdef CUNSEARCH_USE_DOUBLE_PRECISION
	using Real3 = double3;
	#define Real3(x, y, z) make_double3(x, y, z)
	#else
	using Real3 = float3;
	#define Real3(x, y, z) make_float3(x, y, z)
	#endif

	using UInt3 = uint3;
	#define UInt3(x, y, z) make_uint3(x, y, z)
	using Int3 = int3;
	#define Int3(x, y, z) make_int3(x, y, z)

	inline CUDA_METHOD Real3 operator*(Real3 left, Real3 right)
	{
		return Real3(left.x * right.x, left.y * right.y, left.z * right.z);
	}

	inline CUDA_METHOD Real3 operator*(Int3 left, Real3 right)
	{
		return Real3(left.x * right.x, left.y * right.y, left.z * right.z);
	}

	inline CUDA_METHOD Real3 operator*(Real3 left, Int3 right)
	{
		return Real3(left.x * right.x, left.y * right.y, left.z * right.z);
	}

	inline CUDA_METHOD Real3 operator*(Real3 left, Real right)
	{
		return Real3(left.x * right, left.y * right, left.z * right);
	}

	inline CUDA_METHOD Real3 operator*(Real left, Real3 right)
	{
		return Real3(left * right.x, left * right.y, left * right.z);
	}




	inline CUDA_METHOD Real3 operator-(Real3 left, Real3 right)
	{
		return Real3(left.x - right.x, left.y - right.y, left.z - right.z);
	}
	inline CUDA_METHOD Real3 operator+(Real3 left, Real3 right)
	{
		return Real3(left.x + right.x, left.y + right.y, left.z + right.z);
	}

	inline CUDA_METHOD Int3 operator+(Int3 left, Int3 right)
	{
		return Int3(left.x + right.x, left.y + right.y, left.z + right.z);
	}

	inline CUDA_METHOD UInt3 operator+(UInt3 left, UInt3 right)
	{
		return UInt3(left.x + right.x, left.y + right.y, left.z + right.z);
	}


	inline CUDA_METHOD Int3 operator*(Int3 left, int right)
	{
		return Int3(left.x * right, left.y * right, left.z * right);
	}

	inline CUDA_METHOD Int3 operator*(int left, Int3 right)
	{
		return Int3(left * right.x, left * right.y, left * right.z);
	}


	inline CUDA_METHOD void operator-=(Real3 &a, Real3 b)
	{
		a.x -= b.x;
		a.y -= b.y;
		a.z -= b.z;
	}
}

namespace cuNSearch
{
	using Int3 = Int3;
	using UInt3 = UInt3;
	using Real3 = Real3;
}