#ifndef DATAIOH
#define DATAIOH

#include <sys/time.h>
#include <stdio.h>

struct OutputParameters
{
    double slopeRepulsion;
    double forceCutoff;
    double magAttraction;
    double wallStrength;
    double damp;
    double mass;
    int bins;
    int minBinDensity;
    int lowerPressureLimit;
    int upperPressureLimit;
    double optimalCost;
    double nBodyCost;
    double percentDiff;
    int drawn;
    double runTime;
};

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
void logRun(const char filename[MAX_FILE_NAME], OutputParameters params)
{
    FILE *fp = fopen(filename, "a");
    fprintf(fp, "=================================\n");
    fprintf(fp, "\tSlope Repulsion: %f\n", params.slopeRepulsion);
    fprintf(fp, "\tForce Cutoff: %f\n", params.forceCutoff);
    fprintf(fp, "\tMag Attraction: %f\n", params.magAttraction);
    fprintf(fp, "\tWall Strength: %f\n", params.wallStrength);
    fprintf(fp, "\tDamp: %f\n", params.damp);
    fprintf(fp, "\tMass: %f\n\n", params.mass);
    fprintf(fp, "\tBins: %d\n", params.bins);
    fprintf(fp, "\tMin Bin Density: %d\n\n", params.minBinDensity);
    fprintf(fp, "\tLower Pressure Limit: %d\n", params.lowerPressureLimit);
    fprintf(fp, "\tUpper Pressure Limit: %d\n\n", params.upperPressureLimit);
    fprintf(fp, "\tOptimal Cost: %f\n", params.optimalCost);
    fprintf(fp, "\tN-Body Cost: %f\n", params.nBodyCost);
    fprintf(fp, "\tPercent Diff: %f\n", params.percentDiff);
    fprintf(fp, "\tDrawn: %s\n", params.drawn ? "Yes":"No");
    fprintf(fp, "\tRun Time: %.4f s\n", params.runTime/1000000);
    fprintf(fp, "=================================\n");
    fclose(fp);
}


void startTimer(double *timer)
{
	timeval temp;
	gettimeofday(&temp, NULL);
	*timer = (double)(temp.tv_sec * 1000000 + temp.tv_usec);
}

double getTimer(double *timer)
{
	timeval temp;
	gettimeofday(&temp, NULL);
	return((double)(temp.tv_sec * 1000000 + temp.tv_usec) - *timer);
}

void endTimer(double *timer)
{
	timeval temp;
	gettimeofday(&temp, NULL);
	*timer = (double)(temp.tv_sec * 1000000 + temp.tv_usec) - *timer;
}
#endif