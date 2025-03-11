
#include <stdlib.h>
#include <math.h>

#include "utils_timer.h"
#include "utils_serial.h"


int main(int argc, char* argv[])
{
	Timer timer;
	initTimer(&timer);


	int i, j; // grid indexes
	int max_iterations; // number of iterations
	int iteration = 1; // current iteration
	double dt = 100; // largest change in t

	printf("Maximum iterations [100-4000]?\n");
	if (scanf("%d", &max_iterations) != 1) {
		fprintf(stderr, "Error reading input.\n");
		return EXIT_FAILURE;
	}


	startTimer(&timer);


	initialize(); // initialize Temp_last including boundary conditions
	// do until error is minimal or until max steps
	while (dt > MAX_TEMP_ERROR && iteration <= max_iterations) {
		// main calculation: average my four neighbors
		for (i = 1; i <= ROWS; i++) {
			for (j = 1; j <= COLUMNS; j++) {
				Temperature[i][j] = 0.25 * (Temperature_last[i + 1][j] + Temperature_last[i - 1][j] +
					Temperature_last[i][j + 1] + Temperature_last[i][j - 1]);
			}
		}

		dt = 0.0; // reset largest temperature change
		// copy grid to old grid for next iteration and find latest dt
		for (i = 1; i <= ROWS; i++) {
			for (j = 1; j <= COLUMNS; j++) {
				dt = fmax(fabs(Temperature[i][j] - Temperature_last[i][j]), dt);
				Temperature_last[i][j] = Temperature[i][j];
			}
		}
		// periodically print test values
		if ((iteration % 100) == 0) {
			track_progress(iteration);
		}
		iteration++;
	}

	stopTimer(&timer);


	// timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

	double elapsedSec = getElapsedTimeMicroseconds(&timer);


	printf("\nMax error at iteration %d was %f\n", iteration - 1, dt);
	printf("Total time was %f microseconds.\n", elapsedSec);
}
