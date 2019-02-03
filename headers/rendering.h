#ifndef RENDERINGH
#define RENDERINGH

float linearTransform(float x, float shift, float scale, float normal)
{
    return(scale*2.0*(x/normal + shift));
}

void linearTransformPoints(float2 *points, int n, float shift, float scale, float dx, float dy)
{
    for(int i = 0; i < n; i++)
    {
        points[i].x = linearTransform(points[i].x, shift, scale, dx);
        points[i].y = linearTransform(points[i].y, shift, scale, dy);
    }
}

void linearShiftPoints(float2 *points, int n, float shiftX, float shiftY)
{
    for(int i = 0; i < n; i++)
    {
        points[i].x = points[i].x - shiftX;
        points[i].y = points[i].y - shiftY;
    }
}

float linearScalePoints(float2 *points, int n, float scale)
{
    float maxDist = 0.0;
    for(int i = 0; i < n; i++)
    {
        float dist = sqrtf( (points[i].x*points[i].x) + (points[i].y*points[i].y));
        if(dist > maxDist) maxDist = dist;
    }
    
    for(int i = 0; i < n; i++)
    {
        points[i].x = scale*points[i].x/maxDist;
        points[i].y = scale*points[i].y/maxDist;
    }
    
    return(scale/maxDist);
}

#endif