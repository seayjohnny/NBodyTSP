#pragma once
#include <GL/freeglut.h>
#include <GL/gl.h>

#include "structs.cuh"
#include "lintrans.cuh"

#define DIM 768
#define TARGET_FPS 60
#define MAX_STRING 20

RunState rs;
float2 p[N];

double TimePerFrame = 1000/TARGET_FPS;
int CurrentFrame = 0;
int startTime;

int windowId0, windowId1;
int closed = 0;
int clicked = 0;

float scale = 0.9;

float Aged[3] = {1.0, 253.0/255.0, 240.0/255.0};
float LightGray0[3] = {220.0/255.0, 220.0/255.0, 220.0/255.0};
float LightGray1[3] = {206.0/255.0, 206.0/255.0, 206.0/255.0};
float LightGray2[3] = {188.0/255.0, 188.0/255.0, 188.0/255.0};
float DarkGray2[3] = {50.0/255.0, 50.0/255.0, 50.0/255.0};
float DarkGray3[3] = {89.0/255.0, 89.0/255.0, 89.0/255.0};
float Green[3] = {66.0/255.0, 244.0/255.0, 143.0/255.0};
float Red[3] = {244.0/255.0, 65.0/255.0, 65.0/255.0};
float Blue[3] = {109.0/255.0, 158.0/255.0, 235.0/255.0};
float White[3] = {1.0, 1.0, 1.0};
float Black[3] = {0.0, 0.0, 0.0};

void drawGrid(float dx, float dy, float gridLineColor[3], float xAxisColor[3], float yAxisColor[3])
{
    // Draw bounding box
    glLineWidth(0.5);
    (gridLineColor==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(gridLineColor[0], gridLineColor[1], gridLineColor[2]);
    glBegin(GL_LINES);
    for(float x = -scale; x <= scale+0.01; x+=dx)
    {
        glVertex2f(x, scale);
        glVertex2f(x, -scale);
    }

    for(float y = -scale; y <= scale+0.01; y+=dy)
    {
        glVertex2f(-scale, y);
        glVertex2f(scale, y);
    }

    (xAxisColor==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(xAxisColor[0], xAxisColor[1], xAxisColor[2]);
    glVertex2f(0.0, scale);
    glVertex2f(0.0, -scale);

    (yAxisColor==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(yAxisColor[0], yAxisColor[1], yAxisColor[2]);
    glColor3f(0.1, 0.4, 0.4);
    glVertex2f(scale, 0.0);
    glVertex2f(-scale, 0.0);

    glEnd();
    //glFlush();
}

void drawPoint(float2 point, float size, float color[3])
{
    glPointSize(size+10);
    (color==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(color[0], color[1], color[2]);
    glBegin(GL_POINTS);
    glVertex2f(point.x,point.y);
    glEnd();
}

void drawPoints(float2* points, int n, float pointSize, float color[3])
{
    glPointSize(pointSize);
    (color==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(color[0], color[1], color[2]);
    for(int i = 0; i < n; i++)
	{
        glBegin(GL_POINTS);
        glVertex2f(points[i].x,points[i].y);
        glEnd();
    }

}

void drawRect(float2 pos, float2 dim, float shade)
{
    glColor3f(0.31*shade, 0.176*shade, 0.5*shade);
    glBegin(GL_POLYGON);
    glVertex2f(pos.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y+dim.y);
    glVertex2f(pos.x, pos.y+dim.y);
    glEnd();
}

void drawCircle(float2 center, float radius, int segments, float thickness, float color[3])
{
    double pi = 3.14159265358979323846;

    (color==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(color[0], color[1], color[2]);
    glLineWidth(thickness);
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i <= segments; i++)
    {
        float theta = i*2.0*pi/segments;
        glVertex2f(radius*cos(theta)+center.x, radius*sin(theta)+center.y);
    }
    glEnd();
}

void drawText(float x, float y, char string[MAX_STRING])
{
    glPushMatrix();
    glTranslatef(x, y, 0);
    glScalef(1/1250.0, 1/4000.0, 0);
    for( char* c = string; *c; c++)
    {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
    }
    glPopMatrix();
}

void drawDensity(int *density, float2 *densityCenters, float *range, int b, int n)
{
    for(int i = 0; i < b*b; i++)
    {
        if(density[i])
        {
            float xr[2] = {range[i%b], range[i%b+1]};
            float yr[2] = {range[b-1-i/b], range[b-1-i/b+1]};

            drawRect(make_float2(xr[0], yr[0]), make_float2(xr[1]-xr[0], yr[1]-yr[0]), ((float)density[i])/sqrtf(n));

        }
    }
    // linearScalePoints(rs.densityCenters, B*B, scale);
    for(int i = 0; i < b*b; i++)
    {
        if(density[i])
        {
            float color[] = {0.3, 0.3, 1.0}; 
            drawPoint(densityCenters[i], 2.0, color);
        }
    }
}

void drawBubbles(Wall bubbles[B*B])
{
    for(int i = 0; i < B*B;i++)
    {
        if(bubbles[i].strength != 0) drawCircle(bubbles[i].center, bubbles[i].radius*scale, 100, bubbles[i].thickness, Blue);
    }
}

void info()
{
    glutSetWindow(windowId1);
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 1.0, 1.0);
    for(int i = 0;i<N;i++)
    {
        char s[MAX_STRING];
        sprintf(s, "node %d: %f", i, mag(rs.nodes[i].pos));
        float y = -1.8*((float)i/(float)N) + 0.9;
        drawText(-0.9, y, s);
    }

    char pressureText[MAX_STRING];
    sprintf(pressureText, "pressure: %f", rs.pressure);
    drawText(-0.9, -0.9, pressureText);

    glutSwapBuffers();
}

void render()
{
    glutSetWindow(windowId0);
    glClear(GL_COLOR_BUFFER_BIT);
    getNodePositions(p, rs.nodes, N);
    linearScalePoints(p, N, scale);
    linspace(rs.range, -scale, scale, B+1, 1);
    linearScalePoints(rs.densityCenters, B*B, scale);
    drawDensity(rs.densities, rs.densityCenters, rs.range, B, rs.numberOfNodes);


    // drawGrid(dx, dy, LightGray1, LightGray2, LightGray2);
    float2 center = {0,0};
    float innerCircleColor[3];
    float outerCircleColor[3];
    if(rs.innerWall.direction){ for(int i = 0;i<3;i++) innerCircleColor[i] = Green[i];}
    else { for(int i = 0;i<3;i++) innerCircleColor[i] = LightGray2[i]; }

    switch(rs.outerWall.direction)
    {
        case -1:
            for(int i = 0;i<3;i++) outerCircleColor[i] = Red[i];
            break;
        case 0:
            for(int i = 0;i<3;i++) outerCircleColor[i] = LightGray2[i];
            break;
        case 1:
            for(int i = 0;i<3;i++) outerCircleColor[i] = Green[i];
            break;
        default:
            for(int i = 0;i<3;i++) outerCircleColor[i] = Green[i];
            break;
    }
    drawCircle(center, rs.outerWall.radius*scale, 100, rs.outerWall.thickness, outerCircleColor);
    drawCircle(center, rs.innerWall.radius*scale, 100, rs.innerWall.thickness, innerCircleColor);

    for(int i = 0; i < B*B;i++)
    {
        linearScalePoint(&rs.bubbles[i].center, scale);

    }
    drawBubbles(rs.bubbles);
    drawPoints(p, N, 10.0, LightGray0);
    glutSwapBuffers();
}

void drawPath()
{
    glutSetWindow(windowId0);
    glClear(GL_COLOR_BUFFER_BIT);
    getNodeInitPositions(p, rs.nodes, N);
    linearScalePoints(p, N, scale);
    //drawGrid(dx, dy, LightGray1, LightGray2, LightGray2);

    glLineWidth(2.0);
    glColor3f(1.0, 0.3, 0.3);
    glBegin(GL_LINE_LOOP);
    for(int i = 0; i < N;i++)
    {
        glVertex2f(p[rs.nBodyPath[i]].x, p[rs.nBodyPath[i]].y);
    }
    glVertex2f(p[rs.nBodyPath[0]].x, p[rs.nBodyPath[0]].y);
    glEnd();
    drawPoints(p, N, 10.0, LightGray0);
    glutSwapBuffers();
}

void close()
{
    closed = 1;
}

void mouse(int button, int state, int x, int y)
{
    if(button==0)
    {
        clicked = 1;
    }
}

void nBodyRenderInit(int argc, char** argv, RunState runState)
{
    rs = runState;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_MULTISAMPLE);
    glutInitWindowSize(DIM,DIM);
    //glutInitWindowPosition(2800,50);
    glutCreateWindow("Traveling Salesman Problem");
    windowId0 = glutGetWindow();
    glutDisplayFunc(render);
    glutIdleFunc(render);
    glutCloseFunc(close);
    glutMouseFunc(mouse);
    //glClearColor(Aged[0], Aged[1], Aged[2], 1.0);
    glClearColor(0.1725490196, 0.1725490196, 0.1725490196, 1.0);
    // glClearColor(1.0, 1.0, 1.0, 1.0);
    
    glutInitWindowSize(256,768);
    glutInitWindowPosition(glutGet(GLUT_WINDOW_X)+DIM, glutGet(GLUT_WINDOW_Y));
    glutCreateWindow("Info");
    windowId1 = glutGetWindow();
    glutDisplayFunc(info);

    glClearColor(0.1725, 0.1725, 0.1725, 1.0);
    glEnable(GL_MULTISAMPLE_ARB);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
        GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoopEvent();
}

void nBodyRenderUpdate(RunState runState)
{
    rs = runState;
    render();
    glutMainLoopEvent();
    info();
    glutMainLoopEvent();
}

void nBodyDrawPath(RunState runState)
{
    rs = runState;
    drawPath();
    glutMainLoopEvent();
}

bool waitForMouseClick()
{
    glutMainLoopEvent();
    return !clicked;
}

int getWindowState(){ return closed ? 0:1;}
int hasMouseClicked(){ return clicked;}