
#include <stdio.h>


#define COLUMNS 1000
#define ROWS_GLOBAL  1000        // this is a "global" row count
#define NPROC            4        // number of processors
#define ROWS (ROWS_GLOBAL/NPROC)  // number of real local rows

// communication tags
#define DOWN     100
#define UP       101

#define tolerance 0.01

double A_new[ROWS + 2][COLUMNS + 2];
double A[ROWS + 2][COLUMNS + 2];


void initialize(int nproc, int my_rank) {

	double tMin, tMax;  //Local boundary limits
	int i, j;

	for (i = 0; i <= ROWS + 1; i++) {
		for (j = 0; j <= COLUMNS + 1; j++) {
			A[i][j] = 0.0;
		}
	}

	// Local boundry condition endpoints
	tMin = (my_rank) * 100.0 / nproc;
	tMax = (my_rank + 1) * 100.0 / nproc;

	// Left and right boundaries
	for (i = 0; i <= ROWS + 1; i++) {
		A[i][0] = 0.0;
		A[i][COLUMNS + 1] = tMin + ((tMax - tMin) / ROWS) * i;
	}

	// Top boundary (PE 0 only)
	if (my_rank == 0)
		for (j = 0; j <= COLUMNS + 1; j++)
			A[0][j] = 0.0;

	// Bottom boundary (Last PE only)
	if (my_rank == nproc - 1)
		for (j = 0; j <= COLUMNS + 1; j++)
			A[ROWS + 1][j] = (100.0 / COLUMNS) * j;

}


// only called by last PE
void track_progress(int iteration) {
	int i;
	printf("---------- Iteration number: %d ------------\n", iteration);

	// output global coordinates so user doesn't have to understand decomposition
	for (i = 5; i >= 0; i--) {
		printf("[%d,%d]: %5.2f  ", ROWS_GLOBAL - i, COLUMNS - i, A_new[ROWS - i][COLUMNS - i]);
	}
	printf("\n");
}
