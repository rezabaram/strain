#ifndef TIMER_H
#define TIMER_H 
#include<time.h>
class CTimer
	{
	public:
	CTimer():starttime(clock()){
	};
	void reset(){
		starttime=clock();
	}
	double read(){
		return (clock()-starttime)/CLOCKS_PER_SEC;
	}
 	private:
	double starttime;
	};
#endif /* TIMER_H */
