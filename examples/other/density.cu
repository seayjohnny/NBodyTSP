//nvcc SeayJohnnyHW5.cu -o SeayJohnnyHW5 -lglut -lGL -lm

#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "./arrays.h"
#include "./drawing.h"
#include "./density.h"

#define DIM 1024
#define NODES 5000

float2 *nodes = (float2*)malloc((NODES)*sizeof(float2));
__constant__ float2 c_nodes[NODES];

float *pixels;
float *buffer;

int trigger = 1;

float rnd(float x)
{
    return(x*rand() / RAND_MAX);
}
void display()
{
    glClear(GL_COLOR_BUFFER_BIT);

    int b = 32;
    drawDensity(nodes, NODES, b, 1.0);
    drawPoints(nodes, NODES, 5.0);
    drawGrid(2.0/b, 2.0/b, 1.0);
    glFlush();
}

int main(int argc, char** argv)
{ 
    srand( time(NULL) );
    for(int i = 0; i < NODES; i++)
    {
        nodes[i].x = rnd(2.0) - 1.0;
        nodes[i].y = rnd(2.0) - 1.0;
    }
    
    cudaMemcpyToSymbol(c_nodes, nodes, sizeof(float2)*NODES);

    pixels = (float*)malloc(DIM*DIM*3*sizeof(float));
    cudaMalloc(&buffer, DIM*DIM*3*sizeof(float));
    
    // Initialize OpenGL
   	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_MULTISAMPLE);
    glutInitWindowSize(DIM, DIM);

   	glutCreateWindow("GPU | Time to render:\t---");
    glutDisplayFunc(display);

    glClearColor(0.0, 0.0, 0.0, 0.1);
    glEnable(GL_MULTISAMPLE_ARB);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);

    glutMainLoop();

    return(0);

}