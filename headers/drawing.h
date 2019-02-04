#ifndef DRAWINGH
#define DRAWINGH

#include <GL/glut.h>
#include <math.h>

void drawText(float x, float y, char *string)
{
    char *c;
    glRasterPos2f(x, y);
    for(c=string;*c!='\0';c++)
    {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *c);
    }

    glFlush();
}

void drawPoint(float2 point, float size, float color[3])
{
    glPointSize(size+10);
    (color==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(color[0], color[1], color[2]);
    glBegin(GL_POINTS);
    glVertex2f(point.x,point.y);
    glEnd();
    //glFlush();
}

void drawPoints(float2 *points, int n, float size, float color[3])
{
    glPointSize(size+2.5);
	glColor3f(0.0, 0.0, 0.0);
	for(int i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(points[i].x,points[i].y);
		glEnd();
    }

    glPointSize(size);
    (color==NULL) ? glColor3f(1.0, 1.0, 1.0):glColor3f(color[0], color[1], color[2]);
	for(int i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(points[i].x,points[i].y);
		glEnd();
    }
    //glFlush();

}

void drawGrid(float dx, float dy, float scale)
{
    // Draw bounding box
    glLineWidth(1.0);
    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_LINES);
    for(float x = -scale; x <= scale; x+=dx)
    {
        glVertex2f(x, scale);
        glVertex2f(x, -scale);
    }

    for(float y = -scale; y <= scale; y+=dy)
    {
        glVertex2f(-scale, y);
        glVertex2f(scale, y);
    }

    glColor3f(0.1, 0.4, 0.1);
    glVertex2f(0.0, scale);
    glVertex2f(0.0, -scale);

    glColor3f(0.1, 0.4, 0.4);
    glVertex2f(scale, 0.0);
    glVertex2f(-scale, 0.0);

    glEnd();
    //glFlush();
}

void drawRect(float2 pos, float2 dim, float scale)
{
    glColor3f(0.31*scale, 0.176*scale, 0.5*scale);
    glBegin(GL_POLYGON);
    glVertex2f(pos.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y+dim.y);
    glVertex2f(pos.x, pos.y+dim.y);
    glEnd();
    //glFlush();
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
    //glFlush();
}

#endif
