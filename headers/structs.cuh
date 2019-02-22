#pragma once

#define N 48
#include "node.cuh"

struct Wall
{
    double strength;
    int direction;
    float radius;
    float thickness;
} WallDefault = {5000.0, 0, 0, 3.0};

struct RunState
{
    int numberOfNodes;
    Node nodes[N];
    int nBodyPath[N];
    int optimalPath[N];

    float m;
    float p;
    float q;

    float damp;
    float mass;
    
    float pressure;
    float lowerPressureLimit;
    float upperPressureLimit;

    float2 geometricCenter;

    Wall innerWall;
    Wall outerWall;

    float progress;
};