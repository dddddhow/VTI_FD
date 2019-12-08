#include "iostream"
#include "math.h"
#include "stdlib.h"
#include "time.h"
#include "unistd.h"
#include"su_trace_head.h"

#define PI 3.1415926
#define NX 400
#define NZ 400
#define PML 150

#include"model.h"
#include"fd.h"


using namespace std;

int main()
{
    //---------------------------------------------//
    //            计时器--开始                     //
    //---------------------------------------------//
    clock_t start,finish;
    double totaltime;
    start=clock();


    //---------------------------------------------//
    //            正演模拟                         //
    //---------------------------------------------//
	forward();

    //---------------------------------------------//
    //            反演成像                         //
    //---------------------------------------------//



    //---------------------------------------------//
    //            计时器--结束                     //
    //---------------------------------------------//
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    printf("\n此程序的运行时间为: %f 秒\n",totaltime);

    return 0;
}





