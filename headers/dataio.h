#ifndef DATAIOH
#define DATAIOH

#include <stdio.h>

#define MAX_FILE_NAME 20

int getNumberOfLines(const char filename[MAX_FILE_NAME])
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

void loadNodesFromFile(const char filename[MAX_FILE_NAME], float2 *nodes)
{   
    FILE *fp = fopen(filename, "r");
    
    int i = 0;
    while(1)
    {
        fscanf(fp,"%f %f", &nodes[i].x, &nodes[i].y );
        i += 1;
        if(feof(fp))
        {
            break;
        }
    }
    fclose(fp);
}

void loadOptimalPath(const char filename[MAX_FILE_NAME], int *path, int n)
{
    FILE *fp = fopen(filename, "r");
    
    int i = 0;
    while(i < n)
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

#endif