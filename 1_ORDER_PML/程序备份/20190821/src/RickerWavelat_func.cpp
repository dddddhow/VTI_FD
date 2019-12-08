#include"../lib/ricker_wavelet.h"

double ricker_wavelet(int Lw, float freq, float dt, float *signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
	float t, t1;
    for(i=0;i<Lw;i++)
    {
        t=dt*i;
        t1=1/freq;//双边雷克子波
        //t1=0;//单边雷克子波
        signal[i]=(1-2*PI*PI*freq*freq*(t-t1)*(t-t1))*exp(-PI*PI*freq*freq*(t-t1)*(t-t1));
    }
}
