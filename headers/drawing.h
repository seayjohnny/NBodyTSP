#ifndef DRAWINGH
#define DRAWINGH

#include <GL/glut.h>


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

void drawPoints(float2 *points, int n, float size)
{
    
    glPointSize(size+5.0);
	glColor3f(0.0, 0.0, 0.0);
	for(int i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(points[i].x,points[i].y);
		glEnd();
    }

    glPointSize(size);
	glColor3f(1.0, 0.2, 0.2);
	for(int i = 0; i < n; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(points[i].x,points[i].y);
		glEnd();
    }

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
    glFlush();
}

void drawRect(float2 pos, float2 dim, float scale)
{
    glColor3f(scale, scale, 1.5*scale);
    glBegin(GL_POLYGON);
    glVertex2f(pos.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y);
    glVertex2f(pos.x+dim.x, pos.y+dim.y);
    glVertex2f(pos.x, pos.y+dim.y);
    glEnd();
    glFlush();
}

#endif
