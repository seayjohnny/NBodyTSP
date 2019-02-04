#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./headers/dataio.h"
#include "./headers/rendering.h"
#include "./headers/drawing.h"
#include "./headers/density.cuh"

#define DIM 1024
#define BLOCK 1024


#define SLOPE_REPULSION 10000.0
#define FORCE_CUTOFF 10.0
#define MAG_ATRACTION 10.0
#define WALL_STRENGTH 200000.0
#define DAMP 20.0
#define MASS 80.0


#define FP "./datasets/att48/coords.txt"
#define FP_OPT "./datasets/att48/path.txt"
int numberOfNodes = getNumberOfLines(FP);
float2 *orginalCoords = (float2*)malloc((numberOfNodes)*sizeof(float2));
float2 *nodes = (float2*)malloc((numberOfNodes)*sizeof(float2));
float2 *posNbody = (float2*)malloc((numberOfNodes)*sizeof(float2));
float2 *velNbody = (float2*)malloc((numberOfNodes)*sizeof(float2));
float2 *accNbody = (float2*)malloc((numberOfNodes)*sizeof(float2));
float  *massNbody = (float*)malloc((numberOfNodes)*sizeof(float));

int *nBodyPath = (int*)malloc((numberOfNodes)*sizeof(int));
double nBodyCost;

int *optimalPath = (int*)malloc((numberOfNodes)*sizeof(int));
double optimalCost;

double percentDiff;

float4 boundingBox;
float2 geometricCenter;
float normalizingFactor;

int mouse_x, mouse_y;
int displayFlag = 0;

void drawNBodyExtrusion(float2* pos, float innerRadius, float outerRadius, int n);

/* ========================================================================== */
/*                                 DATA SETUP                                 */
/* ========================================================================== */

float2 maxValues(float2 *nodes, int n)
{
    float maxX = nodes[0].x;
    float maxY = nodes[0].y;

    for(int i = 1; i < n; i++)
    {
        if(nodes[i].x > maxX) maxX = nodes[i].x;
        if(nodes[i].y > maxY) maxY = nodes[i].y;
    }
    return(make_float2(maxX, maxY));
}

float2 minValues(float2 *nodes, int n)
{
    float minX = nodes[0].x;
    float minY = nodes[0].y;

    for(int i = 1; i < n; i++)
    {
        if(nodes[i].x < minX) minX = nodes[i].x;
        if(nodes[i].y < minY) minY = nodes[i].y;
    }

    return(make_float2(minX, minY));
}

float2 getGeometricCenter(float2 *nodes, int n)
{
    float2 geometricCenter;
	
	geometricCenter.x = 0.0;
	geometricCenter.y = 0.0;

	for(int i = 0; i < n; i++)
	{
		geometricCenter.x += nodes[i].x;
		geometricCenter.y += nodes[i].y;
	}
	
	geometricCenter.x /= (float)n;
	geometricCenter.y /= (float)n;

	return(geometricCenter);
}

float4 getBoundingBox(float2 *nodes, int n)
{
    float2 min = minValues(nodes, n);
    float2 max = maxValues(nodes, n);

	return(make_float4(min.x, min.y, max.x, max.y));
}

void setNbodyInitialConditions(float2 *node, float2 *pos, float2 *vel, float* mass, int n)
{
	int i;

	for(i = 0; i < n; i++)
	{
		pos[i].x = node[i].x;
		pos[i].y = node[i].y;
		
		vel[i].x = 0.0;
		vel[i].y = 0.0;
		
		mass[i] = MASS;
	}
}

void getNBodyPath(float2 *pos, int* path, int n)
{
	int i;
	double minValue;
	double *angle = (double*)malloc(n*sizeof(double));
    int *used = (int*)malloc(n*sizeof(int));
    double pi = 3.14159265358979323846;
	
	for(i = 0; i < n; i++)
	{
		if(pos[i].x == 0 && pos[i].y == 0)
		{
			angle[i] = 0.0;
		}
		else if(pos[i].x >= 0 && pos[i].y >= 0)
		{
			if(pos[i].x == 0) angle[i] = 90.0;
			else angle[i] = atan(pos[i].y/pos[i].x)*180.0/pi;
		}
		else if(pos[i].x < 0 && pos[i].y >= 0)
		{
			angle[i] = 180.0 - atan(pos[i].y/(-pos[i].x))*180.0/pi;
		}
		else if(pos[i].x <= 0 && pos[i].y < 0)
		{
			if(pos[i].x == 0) angle[i] = 270.0;
			else angle[i] = 180.0 + atan(pos[i].y/pos[i].x)*180.0/pi;
		}
		else
		{
			angle[i] = 360.0 - atan(-pos[i].y/pos[i].x)*180.0/pi;
		}
	}
	
	for(i = 0; i < n; i++)
	{
		used[i] = 0;
	}
	
	for(int k = 0; k < n; k++)
	{
		minValue = 400.0;
		for(i = 0; i < n; i++)
		{
			if(angle[i] < minValue && used[i] == 0)
			{
				minValue = angle[i];
				path[k] = i;
			}
		}
		used[path[k]] = 1;
		//printf("path[%d] = %d\n", k, path[k]);
	}
	
	free(angle);
	free(used);
}

double getPathCost(float2 *node, int *path, int n, int debug)
{
    double cost;
    float dx, dy, temp;

    for(int i = 0; i < n; i++)
    {
        dx = node[path[i%n]].x - node[path[(i+1)%n]].x;
        dy = node[path[i%n]].y - node[path[(i+1)%n]].y;
        temp = sqrt(dx*dx + dy*dy);
        if(debug)
        {
            printf("\t%d (%f, %f) -> %d (%f, %f) = %f\n", path[i%n], node[path[i%n]].x, node[path[i%n]].y, path[(i+1)%n], node[path[(i+1)%n]].x, node[path[(i+1)%n]].y, temp);
            //printf("\t%f\n",temp);

        }
        cost += sqrt(dx*dx + dy*dy);
    }

    return(cost);
}

/* ========================================================================== */
/*                               CUDA FUNCTIONS                               */
/* ========================================================================== */

__global__ void accelerations(float2* nodes, float2* pos, float2* vel, float2* acc, float *mass, float innerRadius, float outerRadius, int n, float dt)
{
    float2 force;
    float2 node, nodePos;
    float dx, dy, dist, edx, edy, edgeLength;
    float radius, forceMag;
    __shared__ float2 shNodes[BLOCK], shPos[BLOCK];

    int id = threadIdx.x;

    force.x = 0.0;
    force.y = 0.0;

    node.x = nodes[id].x;
	node.y = nodes[id].y;
	nodePos.x = pos[id].x;
	nodePos.y = pos[id].y;
		    
    shPos[id] = nodePos;
    shNodes[id] = node;
    __syncthreads();

    for(int i = 0; i < n; i++)	
    {
        if(i != id) 
        {
            dx = shPos[i].x - nodePos.x;
            dy = shPos[i].y - nodePos.y;
            dist = sqrtf(dx*dx + dy*dy);
            
            edx = shNodes[i].x - node.x;
            edy = shNodes[i].y - node.y;
            edgeLength = sqrtf(edx*edx + edy*edy);
            
            if(dist <= edgeLength)
            {
                forceMag = -(edgeLength - dist)*SLOPE_REPULSION;

            }
            else if(edgeLength < dist && dist < FORCE_CUTOFF)
            {
                forceMag =  MAG_ATRACTION/edgeLength;
            }
            else
            {
                forceMag = 0.0;
            }
            
            force.x += forceMag*dx/dist;
            force.y += forceMag*dy/dist;
        }
    }
    
	if(id < n)
	{
		// Forces between node and the walls
		dx = nodePos.x;
		dy = nodePos.y; 
		radius = sqrtf(dx*dx + dy*dy);
	
		if(radius < innerRadius) // Inside inner wall
		{
			forceMag = WALL_STRENGTH*(innerRadius - radius);
			force.x += forceMag*dx/radius;
			force.y += forceMag*dy/radius;
		}
		else if(radius > outerRadius) // Outside outer wall
		{
			forceMag = WALL_STRENGTH*(outerRadius - radius);
			force.x += forceMag*dx/radius;
			force.y += forceMag*dy/radius;
		}
		
		// Adding on damping force.
		force.x += -DAMP*vel[id].x;
		force.y += -DAMP*vel[id].y;
		
		// Creating the accelerations.
	    acc[id].x = force.x/mass[id];
	    acc[id].y = force.y/mass[id];
    }

    __syncthreads();

    if(id < n)
    {
	    vel[id].x += acc[id].x*dt;
		vel[id].y += acc[id].y*dt;
		
		pos[id].x  += vel[id].x*dt;
		pos[id].y  += vel[id].y*dt;
    }

}

double nBodyExtrustionTSP(float2* nodes, float2* pos, float2* vel, float2* acc, float *mass, int n)
{
    float2 *nodesGPU, *posGPU, *velGPU, *accGPU;
    float *massGPU;
    float innerRadius, outerRadius, dr;
    double t;
    float dt = 0.001;
    int drawCount;

    cudaMalloc(&nodesGPU, n * sizeof(float2));
    cudaMalloc(&posGPU, n * sizeof(float2));
    cudaMalloc(&velGPU, n * sizeof(float2));
    cudaMalloc(&accGPU, n * sizeof(float2));
    cudaMalloc(&massGPU, n * sizeof(float));

    cudaMemcpy( nodesGPU, nodes, n * sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy( posGPU, pos, n * sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy( velGPU, vel, n * sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy( accGPU, acc, n * sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy( massGPU, mass, n * sizeof(float), cudaMemcpyHostToDevice);

    innerRadius = 0.0;
    outerRadius = getLargestMagnitude(nodes, n);
    dr = outerRadius/1000;
    drawNBodyExtrusion(pos, innerRadius, outerRadius, n);
    int temp;
    scanf("%d", &temp);
    drawCount = 0;
    while(innerRadius < outerRadius)
    {
        //outerRadius += dr;
        innerRadius += dr;
        t = 0.0;
        while(t < 2.0)
        {
            accelerations<<<1, numberOfNodes>>>(nodesGPU, posGPU, velGPU, accGPU, massGPU, innerRadius, outerRadius, n, dt);
            if(drawCount == 20)
            {
                cudaMemcpy(pos, posGPU, n * sizeof(float2), cudaMemcpyDeviceToHost);
                drawNBodyExtrusion(pos, innerRadius, outerRadius, n);
                drawCount = 0;
            }
            drawCount++;
            t += dt;
        }
    }

    cudaFree(nodesGPU); cudaFree(posGPU); cudaFree(velGPU); cudaFree(accGPU); cudaFree(massGPU);
    return(0);

}

/* ========================================================================== */
/*                              OPENGL FUNCTIONS                              */
/* ========================================================================== */

void drawNBodyExtrusion(float2* pos, float innerRadius, float outerRadius, int n)
{
    int lineAmount = 100;
    
    float2 center = make_float2(0.0, 0.0);
    float normFactor = outerRadius;

    outerRadius /= normFactor;
    innerRadius /= normFactor;

    glClear(GL_COLOR_BUFFER_BIT);
    float innerCircleColor[] = {0.88, 0.61, 0.0};
    float outerCircleColor[] = {0.88, 0.20, 0.0};
    drawCircle(center, innerRadius, lineAmount, 2.0, innerCircleColor);
    drawCircle(center, outerRadius, lineAmount, 2.0, outerCircleColor);
    drawPoints(pos, n, 5.0, NULL);
}

void drawNBodyPath(float2 *nodes, int *path, int n)
{
    glClear(GL_COLOR_BUFFER_BIT);

    glLineWidth(2.0);
    glColor3f(1.0, 0.2, 0.2);
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < n; i++)
    {
        glVertex2f(nodes[path[i]].x, nodes[path[i]].y);
    }
    glEnd();
    glFlush();

    drawPoints(nodes, n, 5.0, NULL);
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    if(displayFlag==0)
    {
        displayFlag = 1;
    }
    else{
        int b = 8;
        float scale = 0.755;
        float colorNodes[3] = {1.0, 0.1, 0.1};
        float colorCenter[3] = {0.0, 1.0, 0.2};
        drawDensity(nodes, numberOfNodes, b, scale);
        drawPoints(nodes, numberOfNodes, 2.5, colorNodes);
        drawGrid(2*scale/b, 2*scale/b, scale);
        drawPoint(geometricCenter, 2.5, colorCenter);
        printf("\t☐  Running n-body extrustion.");
        nBodyExtrustionTSP(nodes, posNbody, velNbody, accNbody, massNbody, numberOfNodes);
        printf("\r\t🗹  N-body extrustion completed.\n");
        getNBodyPath(posNbody, nBodyPath, numberOfNodes);
        drawNBodyPath(nodes, nBodyPath, numberOfNodes);
        printf("\t🗹  N-body path obtained.\n");
        nBodyCost = getPathCost(orginalCoords, nBodyPath, numberOfNodes, 0);
        printf("\t🗹  N-body path cost calculated. Cost = %f\n", nBodyCost);
        percentDiff = 100*(nBodyCost - optimalCost)/optimalCost;
        printf("\t🗹  N-body cost percent difference calculated. Percent Difference = %f%%\n", percentDiff);
    }
}

void mouseMotion(int mx, int my)
{
    mouse_x = mx;
    mouse_y = my;
    //glutPostRedisplay();
}

/* ========================================================================== */
/*                                    MAIN                                    */
/* ========================================================================== */

int main(int argc, char** argv)
{
    setbuf(stdout, NULL);
    /* ----------------------- Loading and Setting up Data ---------------------- */
    loadOptimalPath(FP_OPT, optimalPath, numberOfNodes);
    loadNodesFromFile(FP, orginalCoords);
    printf("\t🗹  Loaded node positions and optimal path.\n");

    optimalCost = getPathCost(orginalCoords, optimalPath, numberOfNodes, 0);
    printf("\t🗹  Obtained optimal cost. Cost = %f\n", optimalCost);

    boundingBox = getBoundingBox(orginalCoords, numberOfNodes);
    printf("\t🗹  Got bounding box.\n");

    float dx = boundingBox.z - boundingBox.x;
    float dy = boundingBox.w - boundingBox.y;


    linearTransformPoints(orginalCoords, nodes, numberOfNodes, -0.5, 1.0, dx, dy);
    printf("\t🗹  Transformed node positions.\n");

    geometricCenter = getGeometricCenter(nodes, numberOfNodes);
    printf("\t🗹  Got geometric center.\n");

	linearShiftPoints(nodes, numberOfNodes, geometricCenter.x, geometricCenter.y);
	geometricCenter = getGeometricCenter(nodes, numberOfNodes);
	printf("\t🗹  Aligned geometric center to origin.\n");

    normalizingFactor = linearScalePoints(nodes, numberOfNodes, 1.0);
	printf("\t🗹  Normalized nodes. Normalizing factor = %f\n", normalizingFactor);

    setNbodyInitialConditions(nodes, posNbody, velNbody, massNbody, numberOfNodes);
    printf("\t🗹  Set n-body initial conditions.\n");

    /* ------------------------------ OpenGL Calls ------------------------------ */
	glutInit(&argc,argv);
	glutInitWindowSize(DIM,DIM);
	glutInitWindowPosition(0,1100);
	glutCreateWindow("Traveling Salesman Problem");
    glutDisplayFunc(display);
    glutPassiveMotionFunc(mouseMotion);
    glClearColor(0.1, 0.1, 0.1, 0.1);
    glEnable(GL_MULTISAMPLE_ARB);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glutMainLoop();
    free(nodes); free(posNbody); free(velNbody); free(accNbody); free(massNbody);
    free(nBodyPath);
}