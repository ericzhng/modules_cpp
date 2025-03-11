/****************************************************************
 * Laplace MPI C Version
 *
 * T is initially 0.0
 * Boundaries are as follows
 *
 *                T                      4 sub-grids
 *   0  +-------------------+  0    +-------------------+
 *      |                   |       |                   |
 *      |                   |       |-------------------|
 *      |                   |       |                   |
 *   T  |                   |  T    |-------------------|
 *      |                   |       |                   |
 *      |                   |       |-------------------|
 *      |                   |       |                   |
 *   0  +-------------------+ 100   +-------------------+
 *      0         T       100
 *
 * Each PE only has a local subgrid.
 * Each PE works on a sub grid and then sends
 * its boundaries to neighbors.
 *
 *******************************************************************/

#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "utils_timer.h"
#include "utils_mpi.h"


int main(int argc, char* argv[])
{
	Timer timer;
	initTimer(&timer);


	int i, j;
	int max_iterations = 1;
	int iteration = 1;
	double dA;

	int        nproc;                // number of PEs
	int        my_rank;           // my PE number
	double     dA_global = 100;       // delta t across all PEs
	MPI_Status status;              // status returned by MPI calls

	// the usual MPI startup routines
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	// verify only NPROC PEs are being used
	if (nproc != NPROC) {
		if (my_rank == 0) {
			printf("This code must be run with %d PEs\n", NPROC);
		}
		MPI_Finalize();
		exit(1);
	}

	// PE 0 asks for input
	if (my_rank == 0) {
		printf("Maximum iterations [100-4000]?\n");
		fflush(stdout); // Not always necessary, but can be helpful

		if (scanf("%d", &max_iterations) != 1) {
			fprintf(stderr, "Error reading input.\n");
			return EXIT_FAILURE;
		}
	}

	// bcast max iterations to other PEs
	MPI_Bcast(&max_iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (my_rank == 0)
		startTimer(&timer);

	initialize(nproc, my_rank);

	while (dA_global > tolerance && iteration <= max_iterations) {

		// main calculation: average my four neighbors
		for (i = 1; i <= ROWS; i++) {
			for (j = 1; j <= COLUMNS; j++) {
				A_new[i][j] = 0.25 * (A[i + 1][j] + A[i - 1][j] +
					A[i][j + 1] + A[i][j - 1]);
			}
		}

		// COMMUNICATION PHASE: send ghost rows for next iteration

		// send bottom real row down
		if (my_rank != nproc - 1) {             //unless we are bottom PE
			MPI_Send(&A_new[ROWS][1], COLUMNS, MPI_DOUBLE, my_rank + 1, DOWN, MPI_COMM_WORLD);
		}

		// receive the bottom row from above into our top ghost row
		if (my_rank != 0) {                  //unless we are top PE
			MPI_Recv(&A[0][1], COLUMNS, MPI_DOUBLE, my_rank - 1, DOWN, MPI_COMM_WORLD, &status);
		}

		// send top real row up
		if (my_rank != 0) {                    //unless we are top PE
			MPI_Send(&A_new[1][1], COLUMNS, MPI_DOUBLE, my_rank - 1, UP, MPI_COMM_WORLD);
		}

		// receive the top row from below into our bottom ghost row
		if (my_rank != nproc - 1) {             //unless we are bottom PE
			MPI_Recv(&A[ROWS + 1][1], COLUMNS, MPI_DOUBLE, my_rank + 1, UP, MPI_COMM_WORLD, &status);
		}

		dA = 0.0;

		for (i = 1; i <= ROWS; i++) {
			for (j = 1; j <= COLUMNS; j++) {
				dA = fmax(fabs(A_new[i][j] - A[i][j]), dA);
				A[i][j] = A_new[i][j];
			}
		}

		// find global dA                                                        
		MPI_Reduce(&dA, &dA_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&dA_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// periodically print test values - only for PE in lower corner
		if ((iteration % 100) == 0) {
			if (my_rank == nproc - 1) {
				track_progress(iteration);
			}
		}
		iteration++;
	}

	// Slightly more accurate timing and cleaner output 
	MPI_Barrier(MPI_COMM_WORLD);

	// PE 0 finish timing and output values
	if (my_rank == 0) {
		stopTimer(&timer);

		// timersub(&stop_time, &start_time, &elapsed_time);
		double elapsedSec = getElapsedTimeMicroseconds(&timer);

		printf("\nMax error at iteration %d was %f\n", iteration - 1, dA_global);
		printf("Total time was %.3f microseconds.\n", elapsedSec);

		fflush(stdout); // Not always necessary, but can be helpful
	}

	MPI_Finalize();
}
