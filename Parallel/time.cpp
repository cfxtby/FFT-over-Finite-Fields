#include "time.hpp"


void Time::start() {
		gettimeofday(&startTime, NULL);
}
void Time::stop() {
	gettimeofday(&endTime, NULL);
	long seconds = endTime.tv_sec - startTime.tv_sec;
	long useconds = endTime.tv_usec - startTime.tv_usec;
	duration = seconds + useconds / 1000000.0;
}
double Time::getDuration() { 
	return duration; 
}
void Time::printTime(){
	printf("%5.6f seconds\n", duration); 
}
