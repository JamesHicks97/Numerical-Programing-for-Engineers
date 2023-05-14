/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 833966
 *   Name        : James Hicks
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	
	char* q2_file = NULL;
	char* q4_file = NULL;
	char* q5_file = NULL;
	double xo;
	char* q6_file = NULL;

	/* TODO: Add timing for each task and output running time in ms */
    
	
	q2_file = argv[1];
	q4_file = argv[2];
	q5_file = argv[3];
	xo = atof(argv[4]);
	q6_file = argv[5];
	
	struct timeval start;
	struct timeval stop;
	
	gettimeofday(&start, NULL);
	shockwave(q2_file);
	gettimeofday(&stop, NULL);
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("TASK 1:  %.2f milliseconds\n", elapsed_ms);
	
	
	
	gettimeofday(&start, NULL);
	linalgbsys(q4_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("TASK 2:  %.2f milliseconds\n", elapsed_ms);
	
	
	gettimeofday(&start, NULL);
	interp(q5_file,xo);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("TASK 3:  %.2f milliseconds\n", elapsed_ms);
	
	
	gettimeofday(&start, NULL);
	waveeqn(q6_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("TASK 4:  %.2f milliseconds\n", elapsed_ms);
	
    
	return (EXIT_SUCCESS);
}
