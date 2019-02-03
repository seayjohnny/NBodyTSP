#ifndef DENSITYH
#define DENSITYH
#include "./arrays.h"
#include "./drawing.h"

__global__ void getDensity(int *filled, float* range, int b, float2* c_nodes, int n)
{
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;
    int id = x + y*gridDim.x*blockDim.x;

    float xr[2] = {range[blockIdx.x], range[blockIdx.x+1]};
    float yr[2] = {range[b-1-blockIdx.y], range[b-1-blockIdx.y+1]};

    for(int i = 0; i < n; i++)
    {
        if(c_nodes[i].x >= xr[0] && c_nodes[i].x <= xr[1])
        {
            if(c_nodes[i].y >= yr[0] && c_nodes[i].y <= yr[1])
            {
                filled[id] += 1;
            }
        }
    }
}


void drawDensity(float2* nodes, int numberOfNodes, int bins, float scale)
{
    int b = bins;
    float* range = (float*)malloc((b+1)*sizeof(float));
    float* range_GPU;
    cudaMalloc(&range_GPU, (b+1)*sizeof(float));
    
    linspace(range, -scale, scale, b+1, 1);
    cudaMemcpy(range_GPU, range, (b+1)*sizeof(float), cudaMemcpyHostToDevice);

    dim3 dimBlock;
    dim3 dimGrid;

    dimBlock.x = 1;
    dimBlock.y = 1;
    dimBlock.z = 1;

    dimGrid.x = b;
    dimGrid.y = b;
    dimGrid.z = 1;

    int *filled_cpu = (int*)malloc(b*b*sizeof(int));
    int *filled;
    cudaMalloc(&filled, sizeof(int)*b*b);

    float2* c_nodes;
    cudaMalloc(&c_nodes, numberOfNodes*sizeof(float2));
    cudaMemcpy(c_nodes, nodes, numberOfNodes*sizeof(float2), cudaMemcpyHostToDevice);
    
    getDensity<<<dimGrid, dimBlock>>>(filled, range_GPU, b, c_nodes, numberOfNodes);
    cudaMemcpy(filled_cpu, filled, sizeof(int)*b*b, cudaMemcpyDeviceToHost);

    for(int i = 0; i < b*b; i++)
    {
        if(filled_cpu[i])
        {
            float xr[2] = {range[i%b], range[i%b+1]};
            float yr[2] = {range[b-1-i/b], range[b-1-i/b+1]};

            drawRect(make_float2(xr[0], yr[0]), make_float2(xr[1]-xr[0], yr[1]-yr[0]), (numberOfNodes/1000 + 1)*((float)filled_cpu[i])/sqrtf(numberOfNodes));
        }
    }
}

#endif