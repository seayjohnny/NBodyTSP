#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./dataio.h"
#include "./rendering.h"
#include "./drawing.h"
#include "./density.h"

#define DIM 1024

int numberOfNodes = getNumberOfLines("data.txt");
float2 *nodes = (float2*)malloc((numberOfNodes)*sizeof(float2));
float4 boundingBox;

int mouse_x, mouse_y;

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


void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    int b = 8;
    float scale = 0.755;
    drawDensity(nodes, numberOfNodes, b, scale);
    drawPoints(nodes, numberOfNodes, 2.5);
    drawGrid(2*scale/b, 2*scale/b, scale);
}

void mouseMotion(int mx, int my)
{
    mouse_x = mx;
    mouse_y = my;
    //glutPostRedisplay();
}

int main(int argc, char** argv)
{
    loadNodesFromFile("data.txt", nodes);    
    printf("🗹  Loaded node positions\n");

    boundingBox = getBoundingBox(nodes, numberOfNodes);
    printf("🗹  Got bounding box\n");

    float dx = boundingBox.z - boundingBox.x;
    float dy = boundingBox.w - boundingBox.y;

    linearTransformPoints(nodes, nodes, numberOfNodes, -0.5, 0.75, dx, dy);
    printf("🗹  Transformed node positions\n");

    
	glutInit(&argc,argv);
	glutInitWindowSize(DIM,DIM);
	glutInitWindowPosition(0,0);
	glutCreateWindow("Traveling Salesman Problem");
    glutDisplayFunc(display);
    glutPassiveMotionFunc(mouseMotion);
    glClearColor(0.1, 0.1, 0.1, 0.1);
    glEnable(GL_MULTISAMPLE_ARB);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
	glutMainLoop();
}