#pragma once

#include <sys/time.h>
#include <stdio.h>
#include "node.cuh"

#ifndef MAX_FILENAME
#define MAX_FILENAME 20
#endif

int getNumberOfNodes(const char filename[MAX_FILENAME])
{
    FILE *fp = fopen(filename, "r");
    int count = 1;
    for(char c = getc(fp); c != EOF; c = getc(fp))
    {
        if(c == '\n') count += 1;
    }

    fclose(fp);
    return count;
}

void loadNodes(Node *nodes, const char filename[MAX_FILENAME])
{
    FILE *fp = fopen(filename, "r");
    
    int i = 0;
    while(1)
    {
        fscanf(fp,"%f %f", &nodes[i].initpos.x, &nodes[i].initpos.y );
        nodes[i].id = i;
        nodes[i].pos.x = nodes[i].initpos.x;
        nodes[i].pos.y = nodes[i].initpos.y;
        nodes[i].vel.x = 0.0;
        nodes[i].vel.y = 0.0;
        nodes[i].acc.x = 0.0;
        nodes[i].acc.y = 0.0;
        nodes[i].mass = 1.0;

        i += 1;
        if(feof(fp))
        {
            break;
        }
    }
    fclose(fp);
}

void loadPath(int* path, const char filename[MAX_FILENAME])
{
    FILE *fp = fopen(filename, "r");
    
    int i = 0;
    while(1)
    {
        fscanf(fp,"%d", &path[i]);
        path[i]--;
        i += 1;
        if(feof(fp))
        {
            break;
        }
    }
    fclose(fp);
}


