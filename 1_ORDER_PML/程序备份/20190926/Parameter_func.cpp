#include"../lib/Parameter_func.h"


double parameter_func(struct PARAMETER* param)
{
    //**模型参数
    param -> NX=100;//x方向样点数
    param -> NZ=100;//z方向样点数
    param -> PML=50;//PML边界厚度
    param -> Nt=2000;//模拟时间样点数
    param -> dx=10.0;//x方向取样
    param -> dz=10.0;//z方向取样
    param -> dt=0.0004;//时间取样

    //**各向异性参数**Thomson表征形式
    param -> epsilon=0.0;
    param -> delt=0.0;
    param -> yita=(param -> epsilon-param -> delt)*1.0/(1+2.0*param -> delt);

    //**有限差分系数
    param -> C1=1.2340911;
    param -> C2=-1.0664985e-01;
    param -> C3=2.3036367e-02;
    param -> C4=-5.3423856e-03;
    param -> C5=1.0772712e-03;
    param -> C6=-1.6641888e-04;
    param -> C7=1.7021711e-005;
    param -> C8=-8.5234642e-007;//差分系数

    //计算空间大小
    param -> Nx=param -> NX+2*param -> PML;//x方向采样点数(x-z方向各留了PML行用于边界处理)
    param -> Nz=param -> NZ+2*param -> PML;

    //子波参数设置
    param -> Lw=param -> Nt;   //子波长度
    param -> freq=28; //主频

    //**观测系统
    param->nx_location= param->Nx/2;
    param->nz_location=param->Nz/2;
    param -> Ns=50; //炮数
    return 0.0;
}
