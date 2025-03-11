
#include <stdio.h>


// size of plate
#define COLUMNS 1000
#define ROWS 1000
// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

double Temperature[ROWS + 2][COLUMNS + 2]; // temperature grid
double Temperature_last[ROWS + 2][COLUMNS + 2]; // temperature grid from last iteration


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize() {
	int i, j;
	for (i = 0; i <= ROWS + 1; i++) {
		for (j = 0; j <= COLUMNS + 1; j++) {
			Temperature_last[i][j] = 0.0;
		}
	}

	// these boundary conditions never change throughout run
	// set left side to 0 and right to a linear increase
	for (i = 0; i <= ROWS + 1; i++) {
		Temperature_last[i][0] = 0.0;
		Temperature_last[i][COLUMNS + 1] = (100.0 / ROWS) * i;
	}

	// set top to 0 and bottom to linear increase
	for (j = 0; j <= COLUMNS + 1; j++) {
		Temperature_last[0][j] = 0.0;
		Temperature_last[ROWS + 1][j] = (100.0 / COLUMNS) * j;
	}
}

// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {
	int i;
	printf("---------- Iteration number: %d ------------\n", iteration);
	for (i = ROWS - 5; i <= ROWS; i++) {
		printf("[%d,%d]: %5.2f ", i, i, Temperature[i][i]);
	}
	printf("\n");
}
