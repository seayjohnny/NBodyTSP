#ifndef RENDERINGH
#define RENDERINGH

float linearTransform(float x, float shift, float scale, float normal)
{
    return(scale*2.0*(x/normal + shift));
}

void linearTransformPoints(float2 *points, float2 *transformedPoints, int n, float shift, float scale, float dx, float dy)
{
    for(int i = 0; i < n; i++)
    {
        transformedPoints[i].x = linearTransform(points[i].x, shift, scale, dx);
        transformedPoints[i].y = linearTransform(points[i].y, shift, scale, dy);
    }
}


#endif