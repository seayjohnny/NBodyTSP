#include <stdio.h>
#include <stdlib.h>

#include "./headers/vecmath.cuh"
#include "./headers/dataio.cuh"
#include "./headers/lintrans.cuh"
#include "./headers/structs.cuh"
#include "./headers/nbodyrender.cuh"

//TODO : nbodyrender.h
//TODO : arrmath.h

#define BLOCK 256
#define FORCE_CUTOFF 100.0

#define DRAW 1

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
    }
}

float getPathCost(RunState rs, RunState os)
{
    float cost;
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
    int n = getNumberOfNodes("./datasets/att48/coords.txt");
    if(n!=N) return 1;

    for(int r = 0; r < 1;r++)
    {
        RunState o_rs;
        RunState h_rs;
        RunState* d_rs;
        cudaMalloc(&d_rs, sizeof(RunState));

        h_rs.numberOfNodes = n;
        h_rs.m = -0.1;
        h_rs.p = 5.0;
        h_rs.q = 7.0;
        h_rs.damp = 20.0;
        h_rs.mass = 80.0;
        h_rs.lowerPressureLimit = 0.5;
        h_rs.upperPressureLimit = 2.0;
        h_rs.innerWall = WallDefault;
        h_rs.innerWall.direction = 1;
        h_rs.outerWall = WallDefault;
        h_rs.outerWall.radius = 1.0;
        h_rs.outerWall.direction = 0;
        h_rs.progress = 0.0;




        loadNodes(h_rs.nodes, "./datasets/att48/coords.txt");
        loadNodes(o_rs.nodes, "./datasets/att48/coords.txt");
        //loadNodes(h_rs.nodes, "./datasets/rand8/coords.txt");
        //loadNodes(o_rs.nodes, "./datasets/rand8/coords.txt");
        
        shiftNodes(h_rs.nodes, n, getGeometricCenter(h_rs.nodes, n)*(-1.0));
        normalizeNodePositions(h_rs.nodes, n);
        resetNodeInitPositions(h_rs.nodes, n);
        
        if(DRAW)
        {
            nBodyRenderInit(argc, argv, h_rs);
            while(waitForMouseClick());
        }
        
        cudaMemcpy(d_rs, &h_rs, sizeof(RunState), cudaMemcpyHostToDevice);

        float innerRadius = 0.0;
        float outerRadius = 1.0;
        float dt = 0.001;
        float dr = outerRadius/100000;
        float t = 0.0;
        int drawCount = 0;
        while(innerRadius < outerRadius-dr)
        {
            t = 0.0;
            while(t < 1.0)
            {
                adjustForPressure<<<1,1>>>(d_rs);
                nBodyStep<<<1,n>>>(d_rs, dt, dr);
                t += dt;
                if(drawCount == 100 && DRAW)
                {
                    cudaMemcpy(&h_rs, d_rs, sizeof(RunState), cudaMemcpyDeviceToHost); //TODO: Stream this maybe
                    cudaDeviceSynchronize();
                    //normalizeNodePositions(h_rs.nodes, n);
                    //normalizeCircles(&h_rs);
                    if(getWindowState()) nBodyRenderUpdate(h_rs);
                    drawCount = 0;
                }
                if(DRAW) drawCount++;
            }
            cudaMemcpy(&innerRadius, &d_rs->innerWall.radius, sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&outerRadius, &d_rs->outerWall.radius, sizeof(float), cudaMemcpyDeviceToHost);
        }
        getNBodyPath(h_rs.nBodyPath, h_rs.nodes, n);
        float cost = getPathCost(h_rs, o_rs);
        float actualCost = 33523.708;
        float percentDiff = 100*(cost - actualCost)/actualCost;
        printf("N-body path cost: %f\n", cost);
        printf("Optimal path cost: %f\n", actualCost);
        printf("Percent difference: %f%%\n", percentDiff);
        if(DRAW) while(getWindowState()) nBodyDrawPath(h_rs);
    }

    return(0);
}