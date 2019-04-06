#include <stdio.h>
#include <stdlib.h>

#include "./headers/vecmath.cuh"
#include "./headers/dataio.cuh"
#include "./headers/lintrans.cuh"
#include "./headers/structs.cuh"
#include "./headers/arrays.h"
#include "./headers/nbodyrender.cuh"
#include "./headers/runtimer.cuh"

//TODO : nbodyrender.h
//TODO : arrmath.h

#define BLOCK 256
#define FORCE_CUTOFF 100.0

#define DRAW 1

__global__ void initBubbles(RunState* rs)
{
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;
    int id = x + y*gridDim.x*blockDim.x;

    if(rs->densities[id] > 3)
    {
        rs->bubbles[id].center = rs->densityCenters[id];
        rs->bubbles[id].strength = 5000.0;
    }
}

__global__ void getDensity(RunState* rs)
{
    int x = threadIdx.x + blockIdx.x*blockDim.x;
    int y = threadIdx.y + blockIdx.y*blockDim.y;
    int id = x + y*gridDim.x*blockDim.x;

    float xr[2] = {rs->range[blockIdx.x], rs->range[blockIdx.x+1]};
    float yr[2] = {rs->range[B-1-blockIdx.y], rs->range[B-1-blockIdx.y+1]};
    float xBar = 0.0;
    float yBar = 0.0;
    if(id == 0){
        for(int i = 0; i < B*B; i++)
        {
            rs->densities[i] = 0;
        }
    }
    __syncthreads();

    for(int i = 0; i < rs->numberOfNodes; i++)
    {
        if(rs->nodes[i].pos.x >= xr[0] && rs->nodes[i].pos.x <= xr[1])
        {
            if(rs->nodes[i].pos.y >= yr[0] && rs->nodes[i].pos.y <= yr[1])
            {
                rs->densities[id] += 1;
                xBar += rs->nodes[i].pos.x;
                yBar += rs->nodes[i].pos.y;
            }
        }
    }
    if(rs->densities[id])
    {
        rs->densityCenters[id].x = xBar/rs->densities[id];
        rs->densityCenters[id].y = yBar/rs->densities[id];
    }
}

__host__ __device__ float getPressure(RunState* rs)
{
    int i;
    float sum, temp;
    
    sum = 0.0;
    for(i = 0; i<rs->numberOfNodes;i++)
    {
        temp = mag(rs->nodes[i].pos) - rs->outerWall.radius;
        if(0 < temp)
        {
            sum += temp;
        }
    }
    return(sum*(rs->outerWall.strength)/(2.0*PI*rs->outerWall.radius)); //PI is defined in vecmath.cuh
}

__global__ void adjustForPressure(RunState* rs)
{
    rs->pressure = getPressure(rs);
    if(rs->pressure < rs->lowerPressureLimit)
    {
        rs->innerWall.direction = 0;
        rs->outerWall.direction = -1;
    }
    else if(rs->pressure < rs->upperPressureLimit)
    {
        rs->innerWall.direction = 1;
        rs->outerWall.direction = 0;
    } 
    else
    {
        rs->innerWall.direction = 0;
        rs->outerWall.direction = 1;
    }
}

__global__ void nBodyStep(RunState* rs, float dt, float dr)
{
    __shared__ float2 shPos[N], shInitPos[N];
    float d, edgeLength;

    int id = threadIdx.x;
    
    if(id < N)
    {
        float2 pos = shPos[id] = rs->nodes[id].pos;
        float2 initPos = shInitPos[id] = rs->nodes[id].initpos;

        float2 force = {0, 0};
        float forceMag, h, c, radius;

        __syncthreads();
        for(int i = 0; i < N; i++)
        {
            if(i != id)
            {
                d = dist(shPos[i], pos);
                edgeLength = dist(shInitPos[i], initPos);

                h = rs->m*(powf(powf(rs->q/rs->p, 1/(rs->q-rs->p))*edgeLength, rs->p))/(1 - rs->p/rs->q);
                c = powf(edgeLength/d, rs->q-rs->p);
                forceMag = (c-1)*h/powf(d, rs->p);
                
                force += make_float2(forceMag*(shPos[i].x - pos.x)/d, forceMag*(shPos[i].y - pos.y)/d);
            }
        }

        radius = mag(pos);
        if(radius < rs->innerWall.radius)
        {
            forceMag = rs->innerWall.strength*(rs->innerWall.radius - radius);
            force += make_float2(forceMag*pos.x/radius, forceMag*pos.y/radius);
        }
        else if(radius > rs->outerWall.radius)
        {
            forceMag = rs->outerWall.strength*(rs->outerWall.radius - radius);
            force += make_float2(forceMag*pos.x/radius, forceMag*pos.y/radius);
        }

        for(int i = 0; i < B*B;i++)
        {
            if(rs->bubbles[i].strength != 0)
            {
                float d = dist(rs->bubbles[i].center, pos);
                if(d < rs->bubbles[i].radius)
                {
                    forceMag = rs->bubbles[i].strength*(rs->bubbles[i].radius - d);
                    force += make_float2(forceMag*pos.x/d, forceMag*pos.y/d);
                }
            }
        }


        force += make_float2(-rs->damp*rs->nodes[id].vel.x, -rs->damp*rs->nodes[id].vel.y);
        if(force.x > FORCE_CUTOFF) force.x = 0;
        if(force.y > FORCE_CUTOFF) force.y = 0;
        rs->nodes[id].acc = force/rs->nodes[id].mass;
       



        __syncthreads();
        rs->nodes[id].vel += rs->nodes[id].acc*dt;
        rs->nodes[id].pos += rs->nodes[id].vel*dt;
    }

    if(id == 0)
    {
        rs->innerWall.radius += dr*rs->innerWall.direction;
        rs->outerWall.radius += dr*rs->outerWall.direction;
        for(int i = 0; i < B*B;i++)
        {
            float d = mag(rs->bubbles[i].center) + rs->bubbles[i].radius;

            if(rs->bubbles[i].strength != 0 && d <=rs->outerWall.radius) 
            {
                rs->bubbles[i].radius += dr*rs->innerWall.direction;
            }
        }
    }
}

float getPathCost(RunState rs, RunState os)
{
    float cost = 0.0;
    float dx, dy;
    int n = rs.numberOfNodes;

    for(int i = 0; i < n; i++)
    {
        dx = os.nodes[rs.nBodyPath[i%n]].initpos.x - os.nodes[rs.nBodyPath[(i+1)%n]].initpos.x;
        dy = os.nodes[rs.nBodyPath[i%n]].initpos.y - os.nodes[rs.nBodyPath[(i+1)%n]].initpos.y;
        cost += sqrt(dx*dx + dy*dy);
    }

    return cost;
}

void normalizeCircles(RunState* rs)
{
    float factor = rs->outerWall.radius;

    rs->innerWall.radius /= factor;
    rs->outerWall.radius /= factor;
}

int main(int argc, char** argv)
{
    int n = getNumberOfNodes("./datasets/bay29/coords.txt");
    if(n!=N) return 1;

    int run = 0;
    
    for(float mv = -0.0020; mv < -0.0010;mv+=0.0001)
    {
        for(float pv = 2.0; pv < 12.0; pv += 0.25)
        {
            for(float qv = pv + 1.0; qv < pv + 4.0; qv += 0.25)
            {
                RunState o_rs;
                RunState h_rs;
                RunState* d_rs;
                cudaMalloc(&d_rs, sizeof(RunState));
                memset(&h_rs, 0, sizeof(h_rs));
                memset(&o_rs, 0, sizeof(o_rs));

                h_rs.numberOfNodes = n;
                h_rs.m = mv;
                h_rs.p = pv;
                h_rs.q = qv;
                h_rs.damp = 20.0;
                h_rs.mass = 80.0;
                h_rs.lowerPressureLimit = 1.0;
                h_rs.upperPressureLimit = 10.0;
                h_rs.innerWall = WallDefault;
                h_rs.innerWall.direction = 1;
                h_rs.innerWall.strength = 5000.0;
                h_rs.outerWall = WallDefault;
                h_rs.outerWall.radius = 1.0;
                h_rs.outerWall.direction = 0;
                h_rs.outerWall.strength = 5000.0;
                h_rs.progress = 0.0;
                memset(h_rs.nBodyPath, 0, n*sizeof(h_rs.nBodyPath[0]));
                linspace(h_rs.range, -1.0, 1.0, B+1, 1);
        
                for(int i = 0; i < B*B; i++)
                {
                    h_rs.bubbles[i] = WallDefault;
                }
                        
                // loadNodes(h_rs.nodes, "./datasets/att48/coords.txt");
                // loadNodes(o_rs.nodes, "./datasets/att48/coords.txt");
                loadNodes(h_rs.nodes, "./datasets/bay29/coords.txt");
                loadNodes(o_rs.nodes, "./datasets/bay29/coords.txt");
                //loadNodes(h_rs.nodes, "./datasets/rand8/coords.txt");
                //loadNodes(o_rs.nodes, "./datasets/rand8/coords.txt");
                //loadNodes(h_rs.nodes, "./datasets/pres8/coords.txt");
                //loadNodes(o_rs.nodes, "./datasets/pres8/coords.txt");
                
                shiftNodes(h_rs.nodes, n, getGeometricCenter(h_rs.nodes, n)*(-1.0));
                normalizeNodePositions(h_rs.nodes, n);
                resetNodeInitPositions(h_rs.nodes, n);
        
                dim3 grid;
                grid.x = B;
                grid.y = B;
                grid.z = 1;
        
                cudaMemcpy(d_rs, &h_rs, sizeof(RunState), cudaMemcpyHostToDevice);
                //getDensity<<<grid, 1>>>(d_rs);
                //initBubbles<<<grid, 1>>>(d_rs);
                //cudaMemcpy(&h_rs, d_rs, sizeof(RunState), cudaMemcpyDeviceToHost);
                if(DRAW)
                {
                    nBodyRenderInit(argc, argv, h_rs);
                    while(waitForMouseClick());
                }

                float innerRadius = 0.0;
                float outerRadius = 1.0;
                float dt = 0.001;
                float dr = outerRadius/100000;
                float t = 0.0;
                int drawCount = 0;
                double timer;
                startTimer(&timer);
                while(innerRadius < outerRadius-dr)
                {
                    t = 0.0;
                    while(t < 1.0)
                    {
                        printf("\r Run %d: Running n-body extrustion. ☐  %.0f%%", run, floor(100*innerRadius/outerRadius));
                        //getDensity<<<grid, 1>>>(d_rs);
                        //adjustForPressure<<<1,1>>>(d_rs);
                        nBodyStep<<<1,n>>>(d_rs, dt, dr);
                        t += dt;
                        if(drawCount == 100 && DRAW)
                        {
                            cudaMemcpy(&h_rs, d_rs, sizeof(RunState), cudaMemcpyDeviceToHost); //TODO: Stream this maybe
                            cudaDeviceSynchronize();
                            normalizeNodePositions(h_rs.nodes, n);
                            normalizeCircles(&h_rs);
                            if(getWindowState()) nBodyRenderUpdate(h_rs);
                            drawCount = 0;
                        }
                        if(DRAW) drawCount++;
                    }
                    cudaMemcpy(&innerRadius, &d_rs->innerWall.radius, sizeof(float), cudaMemcpyDeviceToHost);
                    cudaMemcpy(&outerRadius, &d_rs->outerWall.radius, sizeof(float), cudaMemcpyDeviceToHost);
                }
                printf("\r Run %d: N-body extrustion completed. 🗹            \n", run);
                cudaMemcpy(&h_rs, d_rs, sizeof(RunState), cudaMemcpyDeviceToHost);
                getNBodyPath(h_rs.nBodyPath, h_rs.nodes, n);
                h_rs.nbodyCost = getPathCost(h_rs, o_rs);
                h_rs.optimalCost = 33523.708;
                h_rs.percentDiff = 100*(h_rs.nbodyCost - h_rs.optimalCost)/h_rs.optimalCost;
                //printf("N-body path cost: %f\n", h_rs.nbodyCost);
                //printf("Optimal path cost: %f\n", h_rs.optimalCost);
                //printf("Percent difference: %f%%\n", h_rs.percentDiff);
                endTimer(&timer);
                h_rs.runTime = timer;
                h_rs.drawn = DRAW;
                logRun(&h_rs, "./runlog2.txt");

                cudaFree(d_rs);
                h_rs.nbodyCost = 0;
                //printf("N-body path cost: %f\n", h_rs.nbodyCost);
                h_rs.optimalCost = 0;
                h_rs.percentDiff = 0;
                if(DRAW) while(getWindowState()) nBodyDrawPath(h_rs);
                ++run;
            }
        }
    }
    return(0);
}