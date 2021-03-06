#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "./headers/dataio.h"
#include "./headers/rendering.h"
#include "./headers/drawing.h"
#include "./headers/density.cuh"


#define DIM 1024
#define FP "./datasets/att48/coords.txt"
int numberOfNodes = getNumberOfLines(FP);
float2 *nodes = (float2*)malloc((numberOfNodes)*sizeof(float2));

float4 boundingBox;
float2 geometricCenter;

int mouse_x, mouse_y;

void delay(int number_of_seconds) 
{ 
    // Converting time into milli_seconds 
    int milli_seconds = 1000 * number_of_seconds; 
  
    // Stroing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not acheived 
    while (clock() < start_time + milli_seconds)
        ; 
} 

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

float rnd(float x)
{
    return(x*rand() / RAND_MAX);
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    int b = 8;
    float scale = 1.0;
    float colorNodes[3] = {0.1, 0.1, 0.1};
    float colorCenter[3] = {0.0, 1.0, 0.2};
    //drawDensity(nodes, numberOfNodes, b, scale);
    drawPoints(nodes, numberOfNodes, 3.0, colorNodes);
    //drawGrid(2*scale/b, 2*scale/b, scale);
    //drawPoint(geometricCenter, 2.5, colorCenter);
    glutSwapBuffers();
    delay(50);
    for(int i = 0; i < numberOfNodes; i++)
    {
        nodes[i].x += rnd(0.06) - 0.03;
        nodes[i].y += rnd(0.06) - 0.03;
    }
}

void mouseMotion(int mx, int my)
{
    mouse_x = mx;
    mouse_y = my;
    //glutPostRedisplay();
}

int main(int argc, char** argv)
{
    loadNodesFromFile(FP, nodes);    
    printf("🗹  Loaded node positions\n");

    boundingBox = getBoundingBox(nodes, numberOfNodes);
    printf("🗹  Got bounding box\n");

    float dx = boundingBox.z - boundingBox.x;
    float dy = boundingBox.w - boundingBox.y;

    linearTransformPoints(nodes, nodes, numberOfNodes, -0.5, 0.75, dx, dy);
    printf("🗹  Transformed node positions\n");

    geometricCenter = getGeometricCenter(nodes, numberOfNodes);
    printf("🗹  Got geometric center\n");
    
    /*linearShiftPoints(nodes, numberOfNodes, geometricCenter.x, geometricCenter.y);
    geometricCenter = getGeometricCenter(nodes, numberOfNodes);
    printf("🗹  Aligned geometric center to origin\n");*/
    
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE);
	glutInitWindowSize(DIM,DIM);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Traveling Salesman Problem");
    glutDisplayFunc(display);
    glutPassiveMotionFunc(mouseMotion);
    glClearColor(1.0, 1.0, 1.0, 0.1);
    glEnable(GL_MULTISAMPLE_ARB);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glutMainLoop();
    free(nodes);
}