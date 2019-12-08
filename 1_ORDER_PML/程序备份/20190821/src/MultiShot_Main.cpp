#include "../lib/ShengShen_head.h"
#include "../lib/Parameter_func.h"
#include "../lib/Forward_func.h"
#include "../lib/RickerWavelet_func.h"
#include "../lib/IO_func.h"


using namespace std;
clock_t start,finish;

int main()
{
    //---------------------------------------------//
    //            计时器--开始                     //
    //---------------------------------------------//
    start=clock();


    //---------------------------------------------//
    //            正演模拟                         //
    //---------------------------------------------//
    //*各种参数*//
    PARAMETER param;
    parameter_func(&param);
    //**正演**//
    //forward(NX, NZ, Nt, PML, dx, dz, dt);
    forward_func(&param);

    //---------------------------------------------//
    //            反演成像                         //
    //---------------------------------------------//



    //---------------------------------------------//
    //            计时器--结束                     //
    //---------------------------------------------//
    finish = clock();
    double totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    printf("此程序的运行时间为: %f 秒\n",totaltime);
    return 0;
}





