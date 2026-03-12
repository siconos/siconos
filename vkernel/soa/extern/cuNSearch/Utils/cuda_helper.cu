#include "cuda_helper.h"
#include <cuda_runtime.h>

CUDAException::CUDAException(const char *_const_Message) : std::runtime_error(_const_Message)
{

}

CUDAMallocException::CUDAMallocException(const char *_const_Message) : std::runtime_error(_const_Message)
{

}

CUDAMemCopyException::CUDAMemCopyException(const char *_const_Message) : std::runtime_error(_const_Message)
{

}

void CudaHelper::DeviceSynchronize()
{
	cudaError_t cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess)
	{
		auto temp = cudaGetErrorString(cudaStatus);
		throw CUDAException(temp);
	}
}

void CudaHelper::GetThreadBlocks(unsigned int numberOfElements, unsigned int alignment, /*out*/ unsigned int &numberOfThreadBlocks, /*out*/ unsigned int &numberOfThreads)
{
	numberOfThreads = (numberOfElements / alignment) * alignment;
	numberOfThreadBlocks = (numberOfElements / alignment);
	if (numberOfElements % alignment != 0)
	{
		numberOfThreads += alignment;
		numberOfThreadBlocks++;
	}
}

void CudaHelper::MemcpyHostToDevice(void* host, void* device, size_t size)
{
	cudaError_t cudaStatus = cudaMemcpy(device, host, size, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess)
	{
		throw CUDAMemCopyException("cudaMemcpy() failed!");
	}
}

void CudaHelper::MemcpyDeviceToHost(void* device, void* host, size_t size)
{
	cudaError_t cudaStatus = cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess)
	{
		throw CUDAMemCopyException("cudaMemcpy() failed!");
	}
}

void CudaHelper::CheckLastError()
{
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		auto temp = cudaGetErrorString(cudaStatus);
		throw CUDAException(temp);
	}
}