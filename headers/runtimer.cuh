#pragma once

#include <sys/time.h>

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