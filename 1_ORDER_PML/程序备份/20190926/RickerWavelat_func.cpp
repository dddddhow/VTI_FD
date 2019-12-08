#include"../lib/RickerWavelet_func.h"

//double rickerwavelet_func(int Lw, float freq, float dt, float *signal)
double rickerwavelet_func(struct PARAMETER *param, float *signal)
{
    //雷克子波，freq为频率,f(t)为雷克子波
    int i;
	float t, t1;
    for(i=0;i<param->Lw;i++)
    {
        t=param->dt*i;
        t1=1/param->freq;//双边雷克子波
        //t1=0;//单边雷克子波
        signal[i]=(1-2*PI*PI*param->freq*param->freq*(t-t1)*(t-t1))
            *exp(-PI*PI*param->freq*param->freq*(t-t1)*(t-t1));
    }
}
