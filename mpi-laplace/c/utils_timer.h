
#include <stdio.h>
#include <stdint.h>


#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif


typedef struct {
#ifdef _WIN32
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	LARGE_INTEGER frequency;
#else
	struct timeval start;
	struct timeval end;
#endif
} Timer;


// Initialize the timer
void initTimer(Timer* timer) {
#ifdef _WIN32
	QueryPerformanceFrequency(&timer->frequency);
#endif
}

// Start the timer
void startTimer(Timer* timer) {
#ifdef _WIN32
	QueryPerformanceCounter(&timer->start);
#else
	gettimeofday(&timer->start, NULL);
#endif
}

// Stop the timer
void stopTimer(Timer* timer) {
#ifdef _WIN32
	QueryPerformanceCounter(&timer->end);
#else
	gettimeofday(&timer->end, NULL);
#endif
}

// Get elapsed time in microseconds
double getElapsedTimeMicroseconds(Timer* timer) {
#ifdef _WIN32
	return (double)(timer->end.QuadPart - timer->start.QuadPart) * 1000000.0 / timer->frequency.QuadPart;
#else
	return (timer->end.tv_sec - timer->start.tv_sec) * 1000000.0 + (timer->end.tv_usec - timer->start.tv_usec);
#endif
}
