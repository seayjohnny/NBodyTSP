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

__global__ void nBodyStep(RunState* rs, float dt, float innerRadius)
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
        rs->innerWall.radius = innerRadius;
    }
}

int main(int argc, char** argv)
{
    int n = getNumberOfNodes("./datasets/att48/coords.txt");
    if(n!=N) return 1;

    for(int r = 0; r < 1;r++)
    {
        RunState h_rs;
        RunState* d_rs;
        cudaMalloc(&d_rs, sizeof(RunState));

        h_rs.numberOfNodes = n;
        h_rs.m = -0.0015;
        h_rs.p = 5.0;
        h_rs.q = 7.0;
        h_rs.damp = 20.0;
        h_rs.mass = 80.0;
        h_rs.lowerPressureLimit = 500;
        h_rs.upperPressureLimit = 1000;
        h_rs.innerWall = WallDefault;
        h_rs.outerWall = WallDefault;
        h_rs.outerWall.radius = 1.0;
        h_rs.progress = 0.0;




        loadNodes(h_rs.nodes, "./datasets/att48/coords.txt");

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
        float dt = 0.01;
        float dr = outerRadius/1000;
        float t = 0.0;
        int drawCount = 0;
        while(innerRadius < outerRadius-dr)
        {
            innerRadius += dr;
            t = 0.0;
            while(t < 1.0)
            {
                nBodyStep<<<1,n>>>(d_rs, dt, innerRadius);
                t += dt;
                if(drawCount == 100 && DRAW)
                {
                    cudaMemcpy(&h_rs, d_rs, sizeof(RunState), cudaMemcpyDeviceToHost); //TODO: Stream this maybe
                    cudaDeviceSynchronize();
                    if(getWindowState()) nBodyRenderUpdate(h_rs);
                    drawCount = 0;
                }
                if(DRAW) drawCount++;
            }

        }
        getNBodyPath(h_rs.nBodyPath, h_rs.nodes, n);
        if(DRAW) while(getWindowState()) nBodyDrawPath(h_rs);
    }

    return(0);
}