#ifndef SIGNAL_H
#define SIGNAL_H

#include <signal.h>

extern "C" {
void SignalControlReport(int sig);
}

extern "C" {
	int bSignal=-1;
	void SignalControlReport(int sig) {
		bSignal=-1;
		if(sig==30)bSignal=0;
		if(sig==31)bSignal=1;
		if(sig==32)bSignal=2;
	}
}


#endif /* SIGNAL_H */
