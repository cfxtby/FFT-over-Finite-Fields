#include <sys/time.h>
#include <stdio.h>

class Time {
	private:
		timeval startTime, endTime;
	public:
		double duration;
		void start();
		void stop();
		double getDuration();
		void printTime();
};

